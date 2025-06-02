# setwd("./multip2")
devtools::load_all()
library(rstan)
options(mc.cores = 10)
rstan::rstan_options(auto_write = TRUE)
source("./inst/testing_code/helper.R")


######simulation test##########

# NETWORK INFO: 
n = 20
t = 2
H = t*(t - 1)/2



# true mean of the priors
Sigma = diag(2*t)
mu_0 = rep(0,t)
rho_0 = rep(0,t)
cross_mu_0 = rep(0,H)
cross_rho_0 = rep(0,H)
eta = 2
Sigma_Dist <- "LJK" 
ig_shape = 3
ig_scale = 50



#######extension to L networks ######
L = 20
n_l_seq <- rpois(L, n)
#n_l_seq <- rep(n,L)

#generate a list of L biplex networks: 
L_networks <- list()
for (i in 1:L) {
    params <- list(
    mu = mu_0,
    rho = rho_0,
    cross_mu = cross_mu_0,
    cross_rho = cross_rho_0,
    Sigma = Sigma,
    C =  mvtnorm::rmvnorm(n=n_l_seq[i],mean=rep(0, 2*t),sigma=Sigma)
    )
    dep_net <- simulate_network(n_l_seq[i], t, params)
    names(dep_net) <- paste0("network", 1:2)
    L_networks[[i]] <- dep_net
}

layer_covar = FALSE
group_covar = TRUE
group_covariates <- matrix(runif(2*L), nrow = L, ncol = 2)
D_group_within <- matrix(c(2,1,1,0), nrow = 2, ncol = t)
D_group_cross <- matrix(c(1,0), nrow = 2, ncol = H)
mu_group_covariates_idx <- c(1, 2,1)
rho_group_covariates_idx <- c(2)
cross_mu_group_covariates_idx <- c(1)
cross_rho_group_covariates_idx <- numeric(0)

stan_data_list <- generate_stan_data_list(L, n_l_seq, L_networks, layer_covar = layer_covar)
stan_data <- generate_multilevel_stan_data(stan_data_list)
stan_data <- add_layer_covar(stan_data, stan_data_list)
if (group_covar) {
    stan_data <- add_group_covar(stan_data, group_covariates, D_group_within, D_group_cross,
                              mu_group_covariates_idx, rho_group_covariates_idx,
                              cross_mu_group_covariates_idx, cross_rho_group_covariates_idx)
} else {
    stan_data <- add_group_covar(stan_data)
}

res <- rstan::stan(file = "./inst/stan/multilevel_multiplex_p2_uncentered.stan", data = stan_data, chains = 4, iter = 1000)

saveRDS(res, file = "/home/annihong/projects/multiplex-social-networks/meeting_notes/res_identifiability.rds")

s <-  summary(res)$summary
s[grep("PS", rownames(s)),c(1,4,8,9,10)]



# res <- rstan::stan(file = "/home/annihong/projects/multip2/inst/stan/dev_model.stan", data = stan_data, chains = 1, iter = 10)
# s_old[grep("PS", rownames(s_old)),c(1,4,8,9,10)]
# stan_data$D_within
# stan_data$D_cross






