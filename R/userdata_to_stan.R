#' Helper method to create covariates for a multiplex p2 model
#'
#' This function creates a list of covariates to be used in formate_covariates function. It assumes that all the same covariates are used for each layer of the multiplex network. 
#'
#' @param t An integer indicating the number of layers of networks.
#' @param H An integer indicating the number of pairs of layers in the multiplex network.
#' @param covars A list containing the names of the covariates to include.
#' @param outcome A character string indicating the names of the layers of the multiplex network.
#' @param network_data A list containing the network data.
#' @param actor_data A list containing the actor data.
#' @return A list of covariates of the form covariates = list(is_within = {T, F}, t_or_h, var_name, value) to be used in the format_covariates function.
covariates_helper <- function(t, H, covars, outcome, network_data, actor_data) {
    create_covar <- function(var_name, t_or_h, is_within, value) {
        return(list(var_name=var_name, t_or_h=t_or_h, is_within=is_within, value=value))
    }
    covariates = list()
    pair_names <- get_pair_names(outcome)
    var_name = "sender"
    for (covar in covars[[var_name]]) {
        res = lapply(1:t, create_covar, var_name=var_name, is_within=TRUE, value=actor_data[covar])
        names(res) = paste(var_name, covar, outcome, sep="_")
        covariates = append(covariates, res)
    }

    var_name = "receiver"
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
#' @param covariates list of covariates in the intermediate form: 
#' covariate <- list(is_within = {T, F}, t_or_h, var_name, value)
#'
#' @return A list of covariates formatted to be estimated via multiplex_p2.stan 
format_covariates <- function(t, H, n, covariates) {
    D_within = matrix(rep(0, t*4), ncol=4, dimnames = list(1:t, c("density", "reciprocity", "sender", "receiver"))) # no covariates
    if (H == 0) { #uniplex network
        D_cross = matrix(rep(0, H*2), ncol = 2)
    } else {
        D_cross = matrix(rep(0, H*2), ncol=2, dimnames = list(1:H, c("cross_density", "cross_reciprocity")))
    }

    covariates_stan <- list(
                            sender = array(dim = c(n, 0)),
                            receiver = array(dim = c(n, 0)),
                            density = array(dim = c(n,n, 0)), 
                            reciprocity = array(dim = c(n,n, 0)),
                            cross_density = array(dim = c(n, n, 0)), 
                            cross_reciprocity = array(dim = c(n, n, 0))
                            )

    for (covar_name in names(covariates)) {
        covariate = covariates[[covar_name]]
        var_name = covariate$var_name
        value = covariate$value
        names(value) = covar_name

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
                            alpha_covariates=covariates_stan$sender,
                            beta_covariates=covariates_stan$receiver, 
                            mu_covariates=covariates_stan$density, 
                            rho_covariates=covariates_stan$reciprocity,
                            cross_mu_covariates=covariates_stan$cross_density, 
                            cross_rho_covariates=covariates_stan$cross_reciprocity,
                            alpha_covariates_S = array(apply(covariates_stan$sender, 2, stats::sd)),
                            beta_covariates_S = array(apply(covariates_stan$receiver, 2, stats::sd)),
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

create_stan_data <- function(model_obj) {
    stan_data_net = format_network(model_obj$data$dep_net)
    stan_data_covar = format_covariates(model_obj$t, model_obj$H, model_obj$n, model_obj$model$covar)
    stan_data_prior = model_obj$model$prior
    stan_data = c(stan_data_net, stan_data_covar, stan_data_prior)
    return(stan_data)
}

# default_priors <- function(stan_data, outcome, covars) {

#     t = stan_data$T
#     H = stan_data$H
#     n = stan_data$n
#     pairs=unlist(get_pair_names(outcome)) 

#     priors <- list(
#             mu_mean_prior = array(0, t), 
#             mu_sd_prior = array(10, t), 
#             rho_mean_prior = array(0, t),
#             rho_sd_prior = array(10, t),
#             cross_mu_mean_prior = array(0, H),
#             cross_mu_sd_prior = array(10, H),
#             cross_rho_mean_prior = array(0, H),
#             cross_rho_sd_prior = array(10, H),
#             mu_covariates_mean_prior = array(0, length(stan_data$mu_covariates_S)),
#             mu_covariates_sd_prior = 10/stan_data$mu_covariates_S,
#             rho_covariates_mean_prior = rep(0, length(stan_data$rho_covariates_S)),
#             rho_covariates_sd_prior = 10/stan_data$rho_covariates_S,
#             cross_mu_covariates_mean_prior = array(0, length(stan_data$cross_mu_covariates_S)),
#             cross_mu_covariates_sd_prior = 10/stan_data$cross_mu_covariates_S,
#             cross_rho_covariates_mean_prior = array(0, length(stan_data$cross_rho_covariates_S)),
#             cross_rho_covariates_sd_prior = 10/stan_data$cross_rho_covariates_S,
#             alpha_covariates_mean_prior = array(0, length(stan_data$alpha_covariates_S)),
#             alpha_covariates_sd_prior = 10/stan_data$alpha_covariates_S,
#             beta_covariates_mean_prior = array(0, length(stan_data$beta_covariates_S)),
#             beta_covariates_sd_prior = 10/stan_data$beta_covariates_S,
#             scale_alpha_prior = 3,
#             scale_beta_prior = 50,
#             LJK_eta_prior = 2
#         )
        
#     priors[c("mu_mean_prior", "mu_sd_prior", "rho_mean_prior", "rho_sd_prior")] <- lapply(priors[c("mu_mean_prior", "mu_sd_prior", "rho_mean_prior", "rho_sd_prior")], function(x) setNames(x, outcome))
#     priors[c("cross_mu_mean_prior", "cross_mu_sd_prior", "cross_rho_mean_prior", "cross_rho_sd_prior")] <- lapply(priors[c("cross_mu_mean_prior", "cross_mu_sd_prior", "cross_rho_mean_prior", "cross_rho_sd_prior")], function(x) setNames(x, pairs))
#     priors[c("mu_covariates_mean_prior", "mu_covariates_sd_prior")] <- lapply(priors[c("mu_covariates_mean_prior", "mu_covariates_sd_prior")], function(x) setNames(x, dimnames(stan_data[["mu_covariates"]])[[3]]))
#     priors[c("rho_covariates_mean_prior", "rho_covariates_sd_prior")] <- lapply(priors[c("rho_covariates_mean_prior", "rho_covariates_sd_prior")], function(x) setNames(x, dimnames(stan_data[["rho_covariates"]])[[3]]))
#     priors[c("cross_mu_covariates_mean_prior", "cross_mu_covariates_sd_prior")] <- lapply(priors[c("cross_mu_covariates_mean_prior", "cross_mu_covariates_sd_prior")], function(x) setNames(x, dimnames(stan_data[["cross_mu_covariates"]])[[3]]))
#     priors[c("cross_rho_covariates_mean_prior", "cross_rho_covariates_sd_prior")] <- lapply(priors[c("cross_rho_covariates_mean_prior", "cross_rho_covariates_sd_prior")], function(x) setNames(x, dimnames(stan_data[["cross_rho_covariates"]])[[3]]))
#     priors[c("alpha_covariates_mean_prior", "alpha_covariates_sd_prior")] <- lapply(priors[c("alpha_covariates_mean_prior", "alpha_covariates_sd_prior")], function(x) setNames(x, dimnames(stan_data[["alpha_covariates"]])[[2]]))
#     priors[c("beta_covariates_mean_prior", "beta_covariates_sd_prior")] <- lapply(priors[c("beta_covariates_mean_prior", "beta_covariates_sd_prior")], function(x) setNames(x, dimnames(stan_data[["beta_covariates"]])[[2]]))
#     return(priors)
# }

default_prior_empty_mdl <- function(t, H, dep_lab, pair_lab) {
    prior <- list( 
        mu_mean_prior = array(0, t), 
        mu_sd_prior = array(10, t), 
        rho_mean_prior = array(0, t),
        rho_sd_prior = array(10, t),
        cross_mu_mean_prior = array(0, H),
        cross_mu_sd_prior = array(10, H),
        cross_rho_mean_prior = array(0, H),
        cross_rho_sd_prior = array(10, H),
        scale_alpha_prior = 3,
        scale_beta_prior = 50,
        LJK_eta_prior = 2
    )

    prior[c("mu_mean_prior", "mu_sd_prior", "rho_mean_prior", "rho_sd_prior")] <- lapply(prior[c("mu_mean_prior", "mu_sd_prior", "rho_mean_prior", "rho_sd_prior")], function(x) setNames(x, dep_lab))
    prior[c("cross_mu_mean_prior", "cross_mu_sd_prior", "cross_rho_mean_prior", "cross_rho_sd_prior")] <- lapply(prior[c("cross_mu_mean_prior", "cross_mu_sd_prior", "cross_rho_mean_prior", "cross_rho_sd_prior")], function(x) setNames(x, pair_lab))

    covar_prior <- list(mu_covariates_mean_prior = array(0, 0),
        mu_covariates_sd_prior = array(0, 0),
        rho_covariates_mean_prior = array(0, 0),
        rho_covariates_sd_prior = array(0, 0),
        cross_mu_covariates_mean_prior = array(0, 0),
        cross_mu_covariates_sd_prior = array(0, 0),
        cross_rho_covariates_mean_prior = array(0, 0),
        cross_rho_covariates_sd_prior = array(0, 0),
        alpha_covariates_mean_prior = array(0, 0),
        alpha_covariates_sd_prior = array(0, 0),
        beta_covariates_mean_prior = array(0, 0),
        beta_covariates_sd_prior = array(0, 0))
    return(append(prior, covar_prior))
}

