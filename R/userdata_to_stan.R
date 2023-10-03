covariates_helper <- function(t, H, covars, outcome, network_data, actor_data) {
    #govt_s = list(var_name = "senders", t_or_h = 1, is_within = T, value = types)
    create_covar <- function(var_name, t_or_h, is_within, value) {
        return(list(var_name=var_name, t_or_h=t_or_h, is_within=is_within, value=value))
    }
    covariates = list()
    pair_names <- get_pair_names(outcome)
    var_name = "senders"
    for (covar in covars[[var_name]]) {
        res = lapply(1:t, create_covar, var_name=var_name, is_within=TRUE, value=actor_data[covar])
        names(res) = paste(var_name, covar, outcome, sep="_")
        covariates = append(covariates, res)
    }

    var_name = "receivers"
    for (covar in covars[[var_name]]) {
        res = lapply(1:t, create_covar, var_name=var_name, is_within=TRUE, value=actor_data[covar])
        names(res) = paste(var_name, covar, outcome, sep="_")
        covariates = append(covariates, res)
    }

    var_name = "density"
    for (covar in covars[[var_name]]) {
        res = lapply(1:t, create_covar, var_name=var_name, is_within=TRUE, value=network_data[[covar]])
        names(res) = paste(var_name, covar, outcome, sep="_")
        covariates = append(covariates, res)
    }

    var_name = "reciprocity"
    for (covar in covars[[var_name]]) {
        res = lapply(1:t, create_covar, var_name=var_name, is_within=TRUE, value=network_data[[covar]])
        names(res) = paste(var_name, covar, outcome, sep="_")
        covariates = append(covariates, res)
    }

    var_name = "cross_density"
    for (covar in covars[[var_name]]) {
        res = lapply(1:H, create_covar, var_name=var_name, is_within=FALSE, value=network_data[[covar]])
        names(res) = paste(var_name, covar, pair_names, sep="_")
        covariates = append(covariates, res)
    }

    var_name = "cross_reciprocity"
    for (covar in covars[[var_name]]) {
        res = lapply(1:H, create_covar, var_name=var_name, is_within=FALSE, value=network_data[[covar]])
        names(res) = paste(var_name, covar, pair_names, sep="_")
        covariates = append(covariates, res)
    }

    return(covariates)
}

#' Format covariates data for stan estimation
#' 
#' @param t T: number of layers in the multiplex network
#' @param H number of pairs of layers
#' @param n number of actors
#' @param covars 
#'
#' @return A list of covariates formatted to be estimated via multiplex_p2.stan 
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

    for (covar_name in names(covariates)) {
        covariate = covariates[[covar_name]]
        var_name = covariate$var_name
        value = covariate$value
        names(value) <- covar_name

        if (covariate$is_within) {
            t = covariate$t_or_h
            D_within[t, var_name] = D_within[t, var_name] + 1
            if (sum(dim(value)) == 2*n) { #dyadic covar
                dim(value) <- c(n, n, 1)
                dimnames(value)[[3]] <- covar_name
                covariates_stan[[var_name]] <- abind::abind(covariates_stan[[var_name]], value, along = 3)
            } else { # actor covar 
                covariates_stan[[var_name]] <- cbind(covariates_stan[[var_name]], value)
            }
            
        } else {
            h = covariate$t_or_h
            D_cross[h, var_name] = D_cross[h, var_name] + 1
            dim(value) <- c(n, n, 1)
            dimnames(value)[[3]] <- covar_name
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
#'
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



