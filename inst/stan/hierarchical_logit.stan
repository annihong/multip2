data {
    int<lower=0> N; //number of observations
    int<lower=1> K; // number of groups
    int<lower=0,upper=1> y[N]; // binary outcomes
    int<lower=1,upper=K> group_idx[N];// group index for each observation 
    int<lower=0, upper=1> sample_prior; // 1 mean sample from prior, 0 mean sample from posterior
}
parameters {
    real mu; //overall intercept 
    real<lower=0> sigma; // standard deviation of the random group-level effects
    real beta[K]; // group-level random effects
}
model {
    // Priors
    mu ~ normal(0, 10);
    sigma ~ inv_gamma(1,1);;
    // random effects
    for (k in 1:K){
        beta[k] ~ normal(mu, sigma);
    }
    // Likelihood
    if (sample_prior != 1){
        for (i in 1:N){
            y[i] ~ bernoulli_logit(mu + beta[group_idx[i]]);
        }
    }  
}
generated quantities {
    int<lower=0,upper=1> y_tilde[N];
    real mu_ps;
    real mean_beta;

    mean_beta = mean(beta);
    mu_ps = mu + mean_beta;

    for (i in 1:N) {
         y_tilde[i] = bernoulli_logit_rng(mu + beta[group_idx[i]]);
    }
}
