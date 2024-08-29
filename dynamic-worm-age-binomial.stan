// Dynamic model for worm burden
// Binomial likelihood for egg count data

functions {
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

  real marginal_egg_neg(real M, real k, vector sens, int Nx){
    vector[Nx] LL;
    LL[1] = neg_binomial_2_lpmf(0 | M, k); //true egg negative
    for(i in 2:Nx) //false egg negatives
      LL[i] = neg_binomial_2_lpmf((i-1) | M, k) + log(1-sens[(i-1)]);
   return log_sum_exp(LL);
  }
  
  real marginal_egg_pos(int eggs, vector epg, real M, real k, real h, int Nx){
    vector[Nx] LL;
    LL[1] = negative_infinity();
    for(i in 2:Nx)
      LL[i] = neg_binomial_2_lpmf(eggs | epg[(i-1)], h) + neg_binomial_2_lpmf((i-1) | M, k);
    return log_sum_exp(LL);
  }
}

data {
  int<lower=1> N_studies;
  int<lower=1> N_classes;
  int<lower=0> N_eggs;
  array[N_eggs] int<lower=0> eggs;
  array[N_eggs] int<lower=1, upper=N_classes> age_class_egg;
  vector[N_classes] age_class_mid; //mid of age class for both eggs and worms
  //array[N_studies, N_classes] int<lower=0> tested; //number of egg count observations per age class
  //array[N_studies, N_classes] int<lower=0> positive; //results of egg diagnostics (0|1) per age class
  array[N_classes] int<lower=0> tested; //number of egg count observations per age class
  array[N_classes] int<lower=0> positive; //results of egg diagnostics (0|1) per age class
  int<lower=500> max_worm; //maximum number of worms to consider
  real <lower=0> mu; //lifespan of adult worm
  real <lower=1> L1; //worm fecundity parameter for egg counts
  real <lower=0,upper=1> gamma; //worm fecundity parameter for egg counts
}

transformed data {
  real<lower=0> b; //diagnostic sensitivity parameter
  array[N_classes] int<lower=0> negative; //inverse results of egg diagnostics (0|1) per age class
  real<lower=0> h; //negative binomial dispersion of eggs
  vector[max_worm] worm_burden;
  vector[max_worm] sens;
  b = 1.7;  //diagnostic sensitivity parameter
  h = 0.4;
 for(x in 1:max_worm){
    sens[x] = x / (b + x);
    worm_burden[x] = x;
  }
  for(j in 1:N_classes)
    negative[j] = tested[j] - positive[j];
}

parameters {
  //array[N_studies] real<lower=0> alpha; //rate at which new worms are acquired
  //array[N_studies] real<lower=0> beta; //rate at which new worms are acquired
  real<lower=0> alpha; //rate at which new worms are acquired
  real<lower=0> beta; //rate at which new worms are acquired
  //real<lower=0> alpha_mu; //hyper parameter
  //real<lower=0> alpha_var; //hyper parameter
  //real<lower=0> beta_mu; //hyper parameter
  //real<lower=0> beta_var; //hyper parameter
  real<lower=0> k_mu; //mean of k
  real<lower=0> k_var; //mean of k
  array[N_classes] real<lower=0> k; //NB dispersion of worms}
}

transformed parameters {
    vector[N_classes] M_age; //Mean worm burden
    vector[N_classes] prev; //true parasite prevalence
    vector[N_classes] S; //population level sensitivity
    vector[N_classes] p; //obervered parasite prevalence
    vector[N_eggs] egg_loglik;
    vector[N_classes] p_loglik;
    vector[max_worm] epg_expected; //Expected epg for different worm burdens
    //expected epg given worm burden
    epg_expected = get_epg(worm_burden, L1, gamma); // vectorised function
    //get M by age
      M_age = get_M(alpha, beta, mu, age_class_mid);
      for(j in 1:N_classes){
          prev[j] = 1 - exp(neg_binomial_2_lpmf(0 | M_age[j],  k[j])); //true prev
      }
    
      //only for egg count data
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

model {

    target += sum(egg_loglik);
    target += sum(p_loglik);
  
  //priors
  //alpha ~ normal(alpha_mu, alpha_var);
  alpha ~ gamma(2, 1); // broad prior around 2
  //alpha_var ~ exponential(2);
  beta ~ gamma(1, 20);
  //beta_mu ~ gamma(1, 20); // vaguely informative prior around 0.05
  //beta_var ~ exponential(2);
  k_mu ~ gamma(3.3, 10);
  k_var ~ exponential(2);
  //for(j in 1:N_studies)
  k ~ normal(k_mu, k_var);
}
