#' Format covariates data for stan estimation
#' 
#' @param stan_data list: network and covariates data for the stan program
#' @param chains numeric: number of MCMC chains, default = 4
#' @param warmup numeric: number of warmup MCMC steps, iter = 2 * warmup, default = 4000
#' @param thin numeric: TODO
#' @param parallel boolean: whether or not to parallelize the chains
#'
#' @return stan_fit object: results of the fit
#' @export

fit_stan <- function(stan_data, chains=4, warmup=4000, thin=1, parallel = TRUE) {
    if (parallel) {
        options(mc.cores = parallel::detectCores())
    }
    rstan::rstan_options(auto_write = TRUE)
    
    p2_fit <- rstan::stan(
                    file="../inst/stan/multiplex_p2.stan", 
                    data = stan_data,
                    chains = chains,
                    iter = warmup*2,
                    warmup = warmup,
                    thin = thin
                    )

    return(p2_fit)
}

