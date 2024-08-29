// Discrete model for probability of CCA by age
// intercept is correlated

data {
  int<lower=1> N_baseline; //number of baseline groups
  int<lower=1> N_fluke; //number of liver fluke exposed group

  array[N_baseline] int<lower=1> tested_baseline; //number tested per age class
  array[N_baseline] int<lower=0> positive_baseline; //number positive per age class
  array[N_baseline] int<lower=0> age_lower_baseline; //mid point of age class
  array[N_baseline] int<lower=0> age_upper_baseline; //mid point of age class
  
  array[N_fluke] int<lower=1> tested_fluke; //number tested per age class
  array[N_fluke] int<lower=0> positive_fluke; //number positive per age class
  array[N_fluke] int<lower=0> age_lower_fluke; //mid point of age class
  array[N_fluke] int<lower=0> age_upper_fluke; //mid point of age class

  int<lower=max(age_upper_fluke)> max_age;
  matrix[max_age, max_age] Dmat; //distance matrix
}

parameters {
    vector[max_age] alpha; // logit probability by year
    real beta_fluke; //coefficient for fluke infection
    real<lower=0> etasq; //correlation params
    real<lower=0> rhosq; //correlation params
}

transformed parameters {
    vector[max_age] prob_baseline; // log probability by year
    vector[max_age] prob_fluke;
    matrix[max_age, max_age] SIGMA_Dmat; //matrix for correlation
    vector[N_baseline] pr_pos_b;    // collect probability for age interval
    vector[N_fluke] pr_pos_f;
    
    for(i in 1:(max_age-1))
      for (j in (i+1):max_age) {
          SIGMA_Dmat[i,j] = etasq*exp(-rhosq*pow(Dmat[i,j],2));
          SIGMA_Dmat[j,i] = SIGMA_Dmat[i,j];
    }
    for(k in 1:max_age)
        SIGMA_Dmat[k,k] = etasq + 0.01;
  
    // Calculate prob for each age value
    for(a in 1:max_age){
      prob_baseline[a] = inv_logit(alpha[a]);
      prob_fluke[a] = inv_logit(alpha[a] + beta_fluke);
    }

    // Age groups are interval censored
    for(i in 1:N_baseline){
      pr_pos_b[i] = 1-prod(1-prob_baseline[age_lower_baseline[i]:age_upper_baseline[i]]);
    }
    for(i in 1:N_fluke){
      pr_pos_f[i] = 1-prod(1-prob_fluke[age_lower_fluke[i]:age_upper_fluke[i]]);
    }
    //print(pr_pos_b);
    //print(pr_pos_f);
}

model {
  // Likelihood for baseline cca (Malaysia)
  for(j in 1:N_baseline){
      target += binomial_lpmf(positive_baseline[j] | tested_baseline[j], pr_pos_b[j]);
    }
        
  // Likelihood for fluke-associated cca (N.E. Thailand)      
  for(j in 1:N_fluke){
      target += binomial_lpmf(positive_fluke[j] | tested_fluke[j], pr_pos_f[j]);
    }
    
  //priors
  //alpha ~ normal(-18, 5);
  beta_fluke ~ normal(2, 5);
  rhosq ~ cauchy(0, 1);
  etasq ~ cauchy(0, 1);
  alpha ~ multi_normal( rep_vector(-18, max_age), SIGMA_Dmat);
}
