
#' Format covariates data for stan estimation
#' 
#' @param outcome vector: names of the layers of networks in the network_data as the outcome multiplex network
#' @param network_data list of p + t matrices of dim n x n : network outcomes and network covariates, p is the number of network covariates and t is the number of layers in the multiplex network outcome. 
#' @param actor_data n x p data frame: p covariates of length of the number of actors (n) 
#' @return a MultiP2Fit object 
#' @export
#'
#' @examples
#' networks = replicate(4, matrix(data=1:9, nrow=3), simplify = FALSE)
#' names(networks) <-  c("network1", "network2", "network_covariate1",      "network_covariate2")
#' MultiP2Fit(outcome = c("network1", "network2"), 
#' network_data = networks, 
#' actor_data = data.frame(actor_attr1=1:3, actor_attr2=1:3))
#constructor
MultiP2Fit <- function(network_data,
                     outcome,
                     actor_data=data.frame(),
                     senders_covar = NULL,
                     receivers_covar = NULL,
                     density_covar = NULL,
                     reciprocity_covar = NULL,
                     cross_density_covar = NULL, 
                     cross_reciprocity_covar = NULL,
                     custom_covars=NULL, 
                     chains = 4, iter = 2000, warmup = floor(iter/2),
                     thin = 1, seed = sample.int(.Machine$integer.max, 1)) {

    stopifnot(is.character(outcome))
    stopifnot(is.list(network_data))
    stopifnot(is.data.frame(actor_data))

    covars = list(senders = senders_covar,
                receivers = receivers_covar,
                density = density_covar,
                reciprocity = reciprocity_covar,
                cross_density = cross_density_covar, 
                cross_reciprocity = cross_reciprocity_covar)
    

    M = network_data[outcome]
    stan_network_data = format_network(M)
    t = stan_network_data$T
    H = stan_network_data$H
    n = stan_network_data$n

    if (is.null(custom_covars)) {
        covariates = covariates_helper(t, H, covars, network_data, actor_data)
    } else {
        covariates = custom_covars
    }
    
    stan_covar = format_covariates(t, H, n, covariates)
    stan_data = append(stan_network_data, stan_covar)

    options(mc.cores = parallel::detectCores())
    rstan::rstan_options(auto_write = TRUE)
    
    fpath <- system.file("stan", "multiplex_p2.stan", package="multiplexP2")
    p2_fit = NULL
    p2_fit <- rstan::stan(
                    file=fpath, 
                    data = stan_data,
                    chains = chains,
                    iter = iter,
                    warmup = warmup,
                    thin = thin,
                    seed = seed
                    )
    newMultiP2Fit = structure(
                            list(stan_fit = p2_fit),
                            network_data=network_data,
                            actor_data=actor_data,
                            outcome=outcome,
                            covariates=covariates,
                            stan_data = stan_data,
                            class = "MultiP2Fit"
                            )

    return(newMultiP2Fit)
}

# networks = replicate(4, matrix(data=1:9, nrow=3), simplify = FALSE)
# names(networks) <-  c("network1", "network2", "network_covariate1",      "network_covariate2")
# actor_data = data.frame(actor_attr1=1:3, actor_attr2=1:3)

# MultiP2Fit(outcome = c("network1", "network2"), 
# network_data = networks, 
# actor_data = actor_data, senders_covar=c("actor_attr1", "actor_attr2"), cross_density_covar = c("network_covariate1","network_covariate2"))

# covars = list(senders = c("actor_attr1", "actor_attr2"),
#             receivers = c("actor_attr1", "actor_attr2"),
#             density = c("network_covariate1","network_covariate2"),
#             reciprocity = c("network_covariate1","network_covariate2"),
#             cross_density = c("network_covariate1","network_covariate2"), 
#             cross_reciprocity = NULL)

#' @export
print.MultiP2Fit <- function(x) {
    cat(attr(f, "outcome"))
}

#' @export
summary.MultiP2Fit <- function(x) {
    pattern = "^mu|^rho|^cross|fixed"
    s <- rstan::summary(x$stan_fit)$summary
    res <- s[grep(pattern,rownames(s)),1:10]
    print(res)
}