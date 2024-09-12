
# NETWORK INFO: 
n = 30
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

sample_prior <- function(i) {
    Sigma <- multip2::sample_Sigma(eta, ig_shape, ig_scale, t)
    sampled_params <- multip2::sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
    A_bar_prior <- colMeans(sampled_params$C[,1 + 2*(1:t - 1)])
    B_bar_prior <- colMeans(sampled_params$C[,2 + 2*(1:t - 1)])
    sampled_params$mu <- sampled_params$mu + A_bar_prior + B_bar_prior 
    res = sampled_params$mu
    print(res)
    return(res)
}


n_sims = 500
n_cores = 50
out <- parallel::mclapply(
        1:n_sims,
        sample_prior,
        mc.cores=n_cores)
saveRDS(out, file="/home/annihong/projects/simres/analysis_res/prior_ps_mu.Rds")

