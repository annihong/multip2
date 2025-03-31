#setwd("./multip2")
devtools::load_all()
library(rstan)
options(mc.cores = 10)
rstan::rstan_options(auto_write = TRUE)

# rstan::rstan_options(auto_write = auto_write)

# n = 5
# networks = replicate(4, matrix(data=sample(c(0,1), n^2, replace=TRUE, prob=c(0.7,0.3)), nrow=n), simplify = FALSE) #self loops are ignored
# names(networks) <-  c(paste0("network", c("1_1", "1_2", "2_1", "2_2")))
# dep_net <- networks[1:2]
# m_empty <- Mp2Model(dep_net)
# stan_data_1 <- create_stan_data(m_empty)$fit_res$stan_data

# dep_net <- networks[3:4]
# m_empty <- Mp2Model(dep_net)
# stan_data_2 <- create_stan_data(m_empty)$fit_res$stan_data


# stan_data <- list()
# stan_data <- stan_data_1
# stan_data$L = 2
# stan_data$n <- c(stan_data_1$n, stan_data_2$n)
# stan_data$N <- c(stan_data_1$N, stan_data_2$N)
# stan_data$y_obs <- c(stan_data_1$y_obs, stan_data_2$y_obs)

# stan_data$network_sim = FALSE
# stan_data$prior_sim = TRUE

# rstan::stan(file = "/home/annihong/projects/multip2/inst/stan/mutilevel_multiplex_p2.stan", data = stan_data, chains = 1, iter = 10)




######simulation test##########

# NETWORK INFO: 
n = 20
t = 2
H = t*(t - 1)/2


# true mean of the priors
#Sigma = diag(2*t)
mu_0 = rep(0,t)
rho_0 = rep(0,t)
cross_mu_0 = rep(0,H)
cross_rho_0 = rep(0,H)
eta = 2
Sigma_Dist <- "LJK" 
ig_shape = 3
ig_scale = 50

# Sigma <- sample_Sigma(eta, ig_shape, ig_scale, t)
# sampled_params <- sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)

Sigma = diag(1, 2*t)
params <- list(
    mu = mu_0,
    rho = rho_0,
    cross_mu = cross_mu_0,
    cross_rho = cross_rho_0,
    Sigma = Sigma,
    C =  mvtnorm::rmvnorm(n=n,mean=rep(0, 2*t),sigma=Sigma)
)
network1 <- simulate_network(n, t, params)
network2 <- simulate_network(n, t, params)

dep_net <- network1
m_empty <- Mp2Model(dep_net)
stan_data_1 <- create_stan_data(m_empty)$fit_res$stan_data

dep_net <- network2
m_empty <- Mp2Model(dep_net)
stan_data_2 <- create_stan_data(m_empty)$fit_res$stan_data


stan_data <- list()
stan_data <- stan_data_1
stan_data$L = 2
stan_data$n <- c(stan_data_1$n, stan_data_2$n)
stan_data$N <- c(stan_data_1$N, stan_data_2$N)
stan_data$y_obs <- c(stan_data_1$y_obs, stan_data_2$y_obs)

stan_data$network_sim = FALSE
stan_data$prior_sim = FALSE

res <- rstan::stan(file = "/home/annihong/projects/multip2/inst/stan/mutilevel_multiplex_p2.stan", data = stan_data, chains = 4, iter = 10)

saveRDS(res, file = "/home/annihong/projects/simres/multilevel_p2/sandbox_res.rds")



