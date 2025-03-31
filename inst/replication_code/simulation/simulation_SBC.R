library(SBC)
library(Rmpi)
library(future)


# cl <- parallel::makeCluster(11, type = "MPI")
# library(future)
# plan(cluster, workers = cl)
# print(cl)
# mpi.close.Rslaves()
# mpi.quit()

source("/home/annihong/projects/multiplex-social-networks/stan/coding/helper.R")
source("/home/annihong/projects/multiplex-social-networks/stan/coding/simulations/simulation_helper.R")
source("/home/annihong/projects/multiplex-social-networks/stan/coding/stan_helper.R")

sample_Sigma <- function(eta, ig_shape, ig_scale, t) {
    L_corr = LKJ_decomp(eta,t)
    #sigma = sigma_folded_normal(sigma_0,sigma_sd,t)
    sigma = sigma_inver_gamma(ig_shape, ig_scale, t)
    Sigma = Sigma_LJK(sigma, L_corr)
    return(Sigma)
}


#prior info:
n = 30
t = 2
H = 1
mu_0 = rep(0,t)
rho_0 = rep(0,t)
cross_mu_0 = rep(0,H)
cross_rho_0 = rep(0,H)
eta = 2
Sigma_Dist <- "LJK" 
ig_shape = 3
ig_scale = 25
set.seed(10032024)
n_sim = 10 # Number of SBC iterations to run


# Setup caching of results
cache_dir <- "/home/annihong/projects/simres/SBC/rstan_cache"
if(!dir.exists(cache_dir)) {
  dir.create(cache_dir)
}

use_cmdstanr <- getOption("SBC.vignettes_cmdstanr", FALSE) # Set to false to use rstan instead

if(use_cmdstanr) {
  library(cmdstanr)
} else {
  library(rstan)
  rstan_options(auto_write = TRUE)
}

#options(mc.cores = parallel::detectCores())
options(mc.cores = 2)

# Enabling parallel processing via future
library(future)
plan(multisession)

######## Model set-up #########
rstan_model <- rstan::stan_model("/home/annihong/projects/multip2/inst/stan/multiplex_p2_revert.stan")


####### Generator #########
multip2_generator_single <- function(n, t, eta, ig_shape, ig_scale, mu_0, rho_0, cross_mu_0, cross_rho_0){  # sim_N is the number of data points we are generating
    Sigma <- sample_Sigma(eta, ig_shape, ig_scale, t)
    sampled_params <- sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
    simulated_network <- simulate_network(n, t, sampled_params)

    A_bar_prior <- colMeans(sampled_params$C[,1 + 2*(1:t - 1)])
    B_bar_prior <- colMeans(sampled_params$C[,2 + 2*(1:t - 1)])
    PS_mu_prior <- sampled_params$mu + A_bar_prior + B_bar_prior

    M <- simulated_network #network to fit:
    m_fit <- multip2::Mp2Model(M)
    m_fit <-  multip2::fit(m_fit, fit = FALSE)



  list(
    variables = list(
        mu = sampled_params$mu,
        PS_mu = PS_mu_prior,
        rho = sampled_params$rho,
        'cross_mu[1]' = sampled_params$cross_mu,
        'cross_rho[1]' = sampled_params$cross_rho
    ),
    generated = m_fit$fit_res$stan_data
  )
}

######### generate SBC_dataets ########


multip2_generator <- SBC_generator_function(multip2_generator_single, n=n, t=t, eta=eta, ig_shape=ig_shape, ig_scale=ig_scale, mu_0=mu_0, rho_0=rho_0, cross_mu_0=cross_mu_0, cross_rho_0=cross_rho_0)

multip2_dataset <- generate_datasets(multip2_generator, n_sims = n_sim)

####### defining backend ########
if(use_cmdstanr) {
  multip2_backend <- SBC_backend_cmdstan_sample(
    cmdstan_model, iter_warmup = 1000, iter_sampling = 1000, chains = 2)
} else {
  multip2_backend <- SBC_backend_rstan_sample(
    rstan_model, iter = 20, warmup = 10, chains = 2, thin=1, open_progress=TRUE, include=FALSE, par = c("x_beta", "y_tilde", "C"))  
}

###### compute results #######
results <- compute_SBC(multip2_dataset, multip2_backend, 
                    cache_mode = "results", 
                    cache_location = file.path(cache_dir, "results_test"),
                    keep_fits = TRUE)