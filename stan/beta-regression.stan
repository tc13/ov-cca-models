// Beta regression - CCA clonal mutation timings

data {
  int<lower=1> M; // total mutations
  int<lower=1> N; // total people
  int<lower=1> G; // total genes
  array[M] int<lower=1, upper=N> person; // patient id for each mutation
  array[M] int<lower=1, upper=G> gene; // patient id for each mutation
  vector[M] y; //mutation timings (proportion)
  int<lower=0> nX; //number of explanatory variables
  matrix[M, nX] X; //explanatory variable values
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  vector[G] alpha; // intercept for gene
  //vector[N] eta; // individual effect
  vector[nX] beta; // coefficients
  real<lower=0> kappa; //beta precision
  //real eta_mu; 
  //real<lower=0> eta_sd;
}

// Linear vector of predictors
transformed parameters {
 vector[M] mu;
 vector[M] Z;
 Z = X * beta;
  for(i in 1:M){;
    //mu[i] = inv_logit(alpha[gene[i]] + eta[person[i]] + Z[i]);
    mu[i] = inv_logit(alpha[gene[i]] + Z[i]);
  }
}

// The model to be estimated
model {
  for(i in 1:M)
    target += beta_proportion_lpdf(y[i] | mu[i], kappa);
    
    //priors
    kappa ~ exponential(0.5);
    beta ~ normal(0, 4);
    //eta ~ normal(eta_mu, eta_sd);
    //eta_mu ~ normal(0, 3);
    //eta_sd ~ exponential(1);
    alpha ~ normal(0, 5);
}

