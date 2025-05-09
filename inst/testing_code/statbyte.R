library(rstan)
options(mc.cores = 4)
rstan::rstan_options(auto_write = TRUE)
set.seed(123)

#simulate data
n = 100
K = 10
N = n*K
group_idx = rep(1:K, each = n)
#true parameters
mu = 0
sigma = 1
beta_ks = rnorm(K, mu, sigma)

#simulate data
prob_logit = mu + beta_ks[group_idx]
prob = exp(prob_logit)/(1 + exp(prob_logit))
y = rbinom(N, 1, prob)

#use rstan to fit the model
stan_data = list(
    N = N,
    K = K,
    y = y,
    group_idx = group_idx,
    sample_prior = 0
)

rstan_fit <- rstan::stan(file = "/home/annihong/projects/multip2/inst/stan/hierarchical_logit.stan", data = stan_data, chains = 4, iter = 2000, warmup = 1000, refresh = 100, control = list(adapt_delta = 0.95), pars = c(), include = FALSE, init = "random", seed = 123)



s <- rstan::summary(rstan_fit)$summary
s[1:20,]
s[grep("mu", rownames(s)),c(1, 4, 8, 9, 10)]
s[grep("sigma", rownames(s)),c(1, 4, 8, 9, 10)]



draws <- rstan::extract(rstan_fit, permuted = TRUE)
str(draws)
y_tilde <- draws[["y_tilde"]]
simulated_means <- apply(y_tilde, 1, mean)
# Plot the simulated means
plot(density(simulated_means), main = "Density of Simulated Means", xlab = "Mean", ylab = "Density")
# add a vertical line for the observed mean
abline(v = mean(y), col = "red", lwd = 2)

#prior sim:
stan_data_prior = stan_data
stan_data_prior$sample_prior = 1
rstan_fit_prior <- rstan::stan(file = "/home/annihong/projects/multip2/inst/stan/hierarchical_logit.stan", data = stan_data_prior, chains = 4, iter = 2000, warmup = 1000, refresh = 1, control = list(adapt_delta = 0.95), pars = c(), include = FALSE, init = "random", seed = 1234)

prior_draws <- rstan::extract(rstan_fit_prior)
y_tilde_prior <- prior_draws[["y_tilde"]]
simulated_means <- apply(y_tilde_prior, 1, mean)
# Plot the simulated means
plot(density(simulated_means), main = "Density of Prior Simulated Means", xlab = "Mean", ylab = "Density")
# add a vertical line for the observed mean
abline(v = mean(y), col = "red", lwd = 2)