// Time to CCA cases, Thailand

functions {  
  //probability of fluke infection
  real foi_def_integral(int x, real alpha, real beta){
    real y;
    real y0;
    y = -(alpha*exp(-beta*x)*(beta*x+1))/beta^2;
    y0 = -(alpha*exp(-beta*0)*(beta*0+1))/beta^2;
    return(y-y0);
    }
}

data {
  int<lower=1> N;                           // Total number of cases
  int<lower=1> age_groups;                  // Number of age groups
  array[age_groups] int<lower=0> cases;     // cases per age group
  array[age_groups] int<lower=1> age_lower; // upper bound of age class (years)
  array[age_groups] int<lower=1> age_upper; // lower bound of age class (years)
  int<lower=max(age_upper)> max_age;        // oldest considered age (years)
  real<lower=0> foi_alpha;                  // parameter for foi
  real<lower=0> foi_beta;                   // parameter for foi
  real<lower=0> k;                          // dispersion parameter for prob of infection
}

transformed data {
  // probability of fluke infection
  array[max_age] real<lower=0> foi;
  array[max_age] real<lower=0, upper=1> prob_infection;
  
  // calculate probability of fluke infection from foi parameters
  for(a in 1:max_age){
    foi[a] = foi_def_integral(a, foi_alpha, foi_beta);              // cumulative worm exposure
    prob_infection[a] = 1-exp(neg_binomial_2_lpmf(0 | foi[a], k));  // prob >=1 worm by age
  }
}

parameters {
    //specify gamma parameters as mean (1/rate)
    real<lower=0> induction_mean;
    real<lower=2> induction_shape;
    real<lower=0> latent_mean;
    real<lower=2> latent_shape;
}

transformed parameters {
    vector[max_age] lp_year;          // probability of event in year
    vector[(max_age-1)] latent_lp;       // probability of latent event per year
    vector[(max_age-1)] induction_lp;    // probability of induction event
    vector[age_groups] lp_interval;   // probability of event in age groups
    
    latent_lp[1] = gamma_lcdf(1 | latent_shape, latent_shape/latent_mean);
    induction_lp[1] = gamma_lcdf(1 | induction_shape, induction_shape/induction_mean);
    
    for(a in 2:(max_age-1)){
      latent_lp[a] = log_diff_exp(gamma_lcdf(a | latent_shape, latent_shape/latent_mean), 
                                   gamma_lcdf(a-1 | latent_shape, latent_shape/latent_mean));
                                  
      induction_lp[a] = log_diff_exp(gamma_lcdf(a | induction_shape, induction_shape/induction_mean), 
                                   gamma_lcdf(a-1 | induction_shape, induction_shape/induction_mean));
    }
    
    lp_year[1] = negative_infinity();
    lp_year[2] = negative_infinity();
    //lp_year[2] = latent_lp[1] + induction_lp[1];
    
    //iterate over age of cancer
    for(a in 3:max_age){
        vector[(a-1)] marginal;
        //mutation is indexed by i, the "true age" of mutation
        marginal[1] = negative_infinity();
        //marginal[1] = latent_lp[1] + induction_lp[(a-1)];
        //probability for mutation at age m
        for(m in 2:(a-1)){
            //in fluke affected population                
            //vector[(m-1)] infection;
            //for(z in 1:(m-1)){
            //      infection[z] = log(prob_infection[z]) + 
            //                        induction_lp[(m-z)] + 
            //                        latent_lp[(a-m)];
            marginal[m] = induction_lp[m] + latent_lp[(a-m)];
            //marginal[m] = log_sum_exp(infection);
        }
      
      //sum marginal probabilities over values of a
      lp_year[a] = log_sum_exp(marginal);
    }
  
    // Age groups are interval censored
    for(j in 1:age_groups){
      lp_interval[j] = log1m_exp(sum(log1m_exp(lp_year[age_lower[j]:age_upper[j]])));
    }
}

model {
        
  // Likelihood for fluke-associated cca (N.E. Thailand)      
  for(j in 1:age_groups)
      for(i in 1:cases[j])
        target += lp_interval[j];
    
  //priors
    induction_mean ~ gamma(5, 0.23);
    induction_shape ~ gamma(5, 0.5);
    
    latent_mean ~ gamma(5, 0.18);
    latent_shape ~ gamma(5, 0.5);
}

