// Dynamic model for worm burden (2)
// Age dependent FoI fitted to age categories (Medley & Bundy)
// Single parameters across all studies
// Worm dispersion (k) varies by age category
// Also incorporates {0,1} egg diagnostic data
// Author Thomas Crellen - thomas.crellen@glasgow.ac.uk
// June 2024

functions {
  //custom model functions
   vector get_epg(vector x, real L1, real gamma){
    return L1*x^gamma;
  }
  
  real get_S(real M, real k, vector sens, int max_worm){
    vector[max_worm] X;
    real S;
    for(x in 1:max_worm)
      X[x] = sens[x] * exp(neg_binomial_2_lpmf(x | M, k));
    S = sum(X) / (1-exp(neg_binomial_2_lpmf(0 | M, k)));
    return(S);
  }
  
  vector get_M(real alpha, real beta, real mu, vector ages){
    int Nx;
    Nx = size(ages);
    vector[Nx] M;
    for(a in 1:Nx)
        M[a] = (alpha*exp(-mu*ages[a])/(beta-mu)^2)*(1+exp(ages[a]*(mu-beta))*(ages[a]*(mu-beta)-1));
    return M;
  }
  
  real marginal_expul(int w, real M, real k, real r, int Nx){
    vector[Nx] LL;
    for(i in 1:Nx)
      LL[i] = neg_binomial_2_lpmf((i+w-1) | M, k) + binomial_lpmf(w | (i+w-1), r);
    return log_sum_exp(LL);
  }
  
  int marginal_expul_infer(int w, real M, real k, real r, int Nx){
    vector[Nx] LL;
    vector[Nx] pState;
    real ll_max;
    int true_worm;
    for(i in 1:Nx)
      LL[i] = neg_binomial_2_lpmf((i+w-1) | M, k) + binomial_lpmf(w | (i+w-1), r);
    pState = exp(LL - log_sum_exp(LL));
    ll_max = max(pState);
    for(i in 1:Nx){
      if(pState[i]==ll_max)
         true_worm = (i+w-1);
    }
      return true_worm;
  }
  
  real marginal_egg_neg(real M, real k, vector sens, int Nx){
    vector[Nx] LL;
    LL[1] = neg_binomial_2_lpmf(0 | M, k); //true egg negative
    for(i in 2:Nx) //false egg negatives
      LL[i] = neg_binomial_2_lpmf((i-1) | M, k) + log(1-sens[(i-1)]);
   return log_sum_exp(LL);
  }
  
  int marginal_egg_neg_infer(real M, real k, vector sens, int Nx){
    vector[Nx] LL;
    vector[Nx] pState;
    real ll_max;
    int true_worm;
    LL[1] = neg_binomial_2_lpmf(0 | M, k); //true egg negative
    for(i in 2:Nx) //false egg negatives
      LL[i] = log(1-sens[(i-1)]) + neg_binomial_2_lpmf((i-1) | M, k);
    pState = exp(LL - log_sum_exp(LL));
    ll_max = max(pState);
    for(i in 1:Nx){
      if(pState[i] == ll_max)
         true_worm = (i-1);
    }
    return true_worm;
  }
  
  real marginal_egg_pos(int eggs, vector epg, real M, real k, real h, int Nx){
    vector[Nx] LL;
    LL[1] = negative_infinity();
    for(i in 2:Nx)
      LL[i] = neg_binomial_2_lpmf(eggs | epg[(i-1)], h) + neg_binomial_2_lpmf((i-1) | M, k);
    return log_sum_exp(LL);
  }

  int marginal_egg_pos_infer(int eggs, vector epg, real M, real k, real h, int Nx){
    vector[Nx] LL;
    vector[Nx] pState;
    real ll_max;
    int true_worm;
    LL[1] = negative_infinity();
    for(i in 2:Nx)
      LL[i] = neg_binomial_2_lpmf(eggs | epg[(i-1)], h) + neg_binomial_2_lpmf((i-1) | M, k);
    pState = exp(LL - log_sum_exp(LL));
    ll_max = max(pState);
    for(i in 1:Nx){
      if(pState[i] == ll_max)
        true_worm = (i-1);
      }
    return true_worm;
  }
  
}

data {
  int<lower=1> N_worm_studies;
  int<lower=0> N_egg_studies;
  int<lower=0> N_prev_studies;
  int<lower=1> N_worms; //total individuals with individual worm burdens
  int<lower=0> N_eggs; //total individuals with individual worm burdens
  array[N_worms] int<lower=1, upper=N_worm_studies> study_id_worm; //index individual to study
  array[N_eggs] int<lower=1, upper=(N_egg_studies+N_worm_studies)> study_id_egg; //index individual to study
  array[N_worms] int<lower=0> worms; //recovered worms per person
  array[N_eggs] int<lower=0> eggs; //recovered worms per person
  int<lower=1> N_classes; //number of age groups for worm data
  array[N_worms] int<lower=1, upper=N_classes> age_class_worm; //age class per person
  array[N_eggs] int<lower=1, upper=N_classes> age_class_egg; //age class per person
  vector[N_classes] age_class_mid; //mid of age class for worms
  int<lower=500> max_worm; //maximum number of worms to consider
  array[N_classes] int<lower=0> tested; //number of egg count observations per age class
  array[N_classes] int<lower=0> positive; //results of egg diagnostics (0|1) per age class
}

transformed data {
  int<lower=N_worms> N; //overall N for individual counts
  int<lower=1> N_studies; //total number of studies
  real<lower=0> r; //worm recovery from expulsions
  real<lower=0> b; //diagnostic sensitivity parameter
  real<lower=0> h; //negative binomial dispersion of eggs
  array[N_classes] int<lower=0> negative; //inverse results of egg diagnostics (0|1) per age class
  array[N_worms] int delta_worm;
  vector[max_worm] sens;
  vector[max_worm] worm_burden;
  r = 0.44; //worm recovery from expulsions
  b = 1.7;  //diagnostic sensitivity parameter
  h = 0.4;
  N = N_worms + N_eggs;
  N_studies = N_worm_studies + N_egg_studies + N_prev_studies;
  for(x in 1:max_worm){
    sens[x] = x / (b + x);
    worm_burden[x] = x;
  }
  //number of additional worms to consider for expulsion studies
  for(i in 1:N_worms){
    delta_worm[i] = max_worm-worms[i];
  }
  //inverse number of results for binary diagnostics
  for(j in 1:N_classes)
    negative[j] = tested[j] - positive[j];

}
  
parameters {
  real<lower=0> alpha; //parameter for force of infection (varies by study)
  real<lower=0> beta; //parameter for force of infection (constant)
  real<lower=0> mu; //worm death rate (constant across studies)
  array[N_classes] real<lower=0> k; // k varies by age class
  real<lower=0> k_mu; //hyper parameter for k
  real<lower=0> k_var; //hyper parameter for k
  real<lower=0> L1; //parameter for worm fecundity
  real<lower=0, upper=1> gamma; //parameter for worm fecundity
}

transformed parameters {
  vector[N_worms] worm_loglik;
  vector[N_eggs] egg_loglik;
  vector[max_worm] epg_expected; //Expected epg for different worm burdens
  vector[N_classes] M_age;
  vector[N_classes] prev; //true parasite prevalence
  vector[N_classes] S; //population level sensitivity
  vector[N_classes] p; //obervered parasite prevalence
  vector[N_classes] p_loglik;
  epg_expected = get_epg(worm_burden, L1, gamma); // vectorised function
  M_age = get_M(alpha, beta, mu, age_class_mid);
  
    for(j in 1:N_classes){
          prev[j] = 1 - exp(neg_binomial_2_lpmf(0 | M_age[j],  k[j])); //true prev
    }
    
  
   // log-likelihood for worm observations
   for(i in 1:N_worms){
     if(study_id_worm[i]==1){ //autopsy
        worm_loglik[i] = neg_binomial_2_lpmf(worms[i] | M_age[age_class_worm[i]], k[age_class_worm[i]]);
     }else{ // expulsion studies
        worm_loglik[i] = marginal_expul(worms[i], M_age[age_class_worm[i]], k[age_class_worm[i]], r, delta_worm[i]);
        }
      }
  
  //log-likelihood for egg observations
  for(i in 1:N_eggs){
    if(eggs[i]==0){ //zero eggs observed
      egg_loglik[i] = marginal_egg_neg(M_age[age_class_egg[i]], k[age_class_egg[i]], sens, max_worm);
    }else{ //at least one egg observed
      egg_loglik[i] = marginal_egg_pos(eggs[i], epg_expected, M_age[age_class_egg[i]], k[age_class_egg[i]], h, max_worm);
    }
  }
  
  //observed prev - for prev data only
   p_loglik = rep_vector(0, N_classes);
    for(j in 1:N_classes){ 
        S[j] = get_S(M_age[j], k[j], sens, max_worm);
        p[j] = prev[j] * S[j];
        for(i in 1:positive[j])
          p_loglik[j] += bernoulli_lpmf(1 | p[j]);
        for(i in 1:negative[j])
          p_loglik[j] += bernoulli_lpmf(0 | p[j]); 
      }
      
}

// The model to be estimated
model {
  //Increment log liklihood
  target += sum(worm_loglik);
  target += sum(egg_loglik);
  target += sum(p_loglik);

  //priors
  alpha ~ gamma(2, 1); // broad prior around 2
  beta ~ gamma(1, 20); // vaguely informative prior around 0.05
  mu ~ gamma(4, 40);  // vaguely informative prior around 1/10
  k ~ normal(k_mu, k_var);
  k_mu ~ gamma(3.3, 10);
  k_var ~ exponential(2);
  L1 ~ gamma(710, 10); // parameters for worm fecundity
  gamma ~ beta(15, 2); // parameters for worm fecundity
}

generated quantities {
  vector[(N+N_classes)] log_lik; //store log likelihood for model comparison
  vector[N_worms] wb_infered_worms;
  vector[N_eggs] wb_infered_eggs;
  
  for(i in 1:N_worms){
    log_lik[i] = worm_loglik[i];
    //Infer worm burden
    if(study_id_worm[i]==1){ //autopsy
        wb_infered_worms[i] = worms[i];
     }else{ // expulsion studies
        wb_infered_worms[i] = marginal_expul_infer(worms[i], M_age[age_class_worm[i]], k[age_class_worm[i]], r, delta_worm[i]);
        }
  }
  
  for(i in 1:N_eggs){
    log_lik[(i+N_worms)] = egg_loglik[i];
    //infer worm burden
    if(eggs[i]==0){ //zero eggs observed
      wb_infered_eggs[i] = marginal_egg_neg_infer(M_age[age_class_egg[i]], k[age_class_egg[i]], sens, max_worm);
    }else{ //at least one egg observed
      wb_infered_eggs[i] = marginal_egg_pos_infer(eggs[i], epg_expected, M_age[age_class_egg[i]], k[age_class_egg[i]], h, max_worm);
    }
  }
  
  for(j in 1:N_classes)
    log_lik[(N+j)] = p_loglik[j];
}
