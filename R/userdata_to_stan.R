#' Format covariates data for stan estimation
#' 
#' @param t T: number of layers in the multiplex network
#' @param H number of pairs of layers
#' @param n number of actors
#' @param covariates list of all the formatted covariates: list(is_within = {T, F}, t_or_h, var_name, value) (still need a function for this)
#'
#' @return A list of covariates formatted to be estimated via multiplex_p2.stan 
#' @export
#'
#' @examples
#' format_covariates(2, 1, 30, covariates=list())

format_covariates <- function(t, H, n, covariates) {
    D_within = matrix(rep(0, t*4), ncol=4, dimnames = list(1:t, c("density", "reciprocity", "senders", "receivers"))) # no covariates
    if (H == 0) { #uniplex network
        D_cross = matrix(rep(0, H*2), ncol = 2)
    } else {
        D_cross = matrix(rep(0, H*2), ncol=2, dimnames = list(1:H, c("cross_density", "cross_reciprocity")))
    }

    covariates_stan <- list(
                            senders = array(dim = c(n, 0)),
                            receivers = array(dim = c(n, 0)),
                            density = array(dim = c(n,n, 0)), 
                            reciprocity = array(dim = c(n,n, 0)),
                            cross_density = array(dim = c(n, n, 0)), 
                            cross_reciprocity = array(dim = c(n, n, 0))
                            )

    for (covariate in covariates) {
        var_name = covariate$var_name
        value = covariate$value
        if (covariate$is_within) {
            t = covariate$t_or_h
            D_within[t, var_name] = D_within[t, var_name] + 1
            if (sum(dim(value)) == 2*n) { #dyadic covar
                covariates_stan[[var_name]] <- abind::abind(covariates_stan[[var_name]], value, along = 3)
            } else { # actor covar 
                covariates_stan[[var_name]] <- cbind(covariates_stan[[var_name]], value)
            }
            
        } else {
            h = covariate$t_or_h
            D_cross[h, var_name] = D_cross[h, var_name] + 1
            covariates_stan[[var_name]] <- abind::abind(covariates_stan[[var_name]], value, along = 3)
        }
    }

    covariates_data <- list(D_within = D_within,
                            D_cross = D_cross,
                            alpha_covariates=covariates_stan$senders,
                            beta_covariates=covariates_stan$receivers, 
                            mu_covariates=covariates_stan$density, 
                            rho_covariates=covariates_stan$reciprocity,
                            cross_mu_covariates=covariates_stan$cross_density, 
                            cross_rho_covariates=covariates_stan$cross_reciprocity,
                            alpha_covariates_S = array(apply(covariates_stan$senders, 2, stats::sd)),
                            beta_covariates_S = array(apply(covariates_stan$receivers, 2, stats::sd)),
                            mu_covariates_S = array(apply(covariates_stan$density, 3, stats::sd)),
                            rho_covariates_S = array(apply(covariates_stan$reciprocity, 3, stats::sd)),
                            cross_mu_covariates_S = array(apply(covariates_stan$cross_density, 3, stats::sd)),
                            cross_rho_covariates_S = array(apply(covariates_stan$cross_reciprocity, 3, stats::sd))
                            )

    return(covariates_data)

}

#' Format network data for stan estimation
#' 
#' @param M list of the observed networks in adjacency matrix form
#' @return list of formatted network related data to feed to Stan of the form list(n=n, N=N, T=t, H=H, layer_pairs=pairs,y_obs=y, K=K, N_obs = length(y),obs_idx = 1:N)

format_network <- function(M) {
    n = nrow(M[[1]]) # number of actors
    t = length(M) # number of layers
    H = t*(t - 1)/2 # number of pairs
    N = n*(n - 1)/2 # number of dyads
    K = 2^(2*t) # number of possible outcomes on a dyad
    pairs = get_dyads(t) # pairs of networks 
    y = matrix_to_dyads_nd(M) # observed networks as dyadic outcomes

    network_data = list(n=n, N=N, T=t, H=H, layer_pairs=pairs,y_obs=y, K=K, N_obs = length(y),obs_idx = 1:N)
    return(network_data)
}


#' The main function to convert user network input to Stan data
#' 
#' @param M list of the observed networks in adjacency matrix form
#' @param covariates intermediate form of the covariates: covariate <- list(is_within = {T, F}, t_or_h, var_name, value)
#' @return a list in the format for Stan code

stan_data <- function(M, covariates) {
    network_data <- format_network(M)
    covariates_data <- format_covariates(covariates)
    stan_data <- append(network_data, covariates_data)
    return(stan_data)
}

