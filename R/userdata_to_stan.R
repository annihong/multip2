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
format_covariate <- function(t, H, n, covariates, prior, data) {
    dict_var_name_to_prior_name <- list("density" = c("mu_covariates_mean_prior", "mu_covariates_sd_prior"), 
                                        "reciprocity" = c("rho_covariates_mean_prior", "rho_covariates_sd_prior"),
                                        "sender" = c("alpha_covariates_mean_prior", "alpha_covariates_sd_prior"),
                                        "receiver" = c("beta_covariates_mean_prior", "beta_covariates_sd_prior"),
                                        "cross_density" = c("cross_mu_covariates_mean_prior", "cross_mu_covariates_sd_prior"),
                                        "cross_reciprocity" = c("cross_rho_covariates_mean_prior", "cross_rho_covariates_sd_prior"))

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
        #print(covar_name)
        covariate = covariates[[covar_name]]
        var_name = covariate$var_name # sender, receiver, density, reciprocity, cross_density, cross_reciprocity
        covar_lab = covariate$covar_lab
        if (covar_lab %in% names(data$dyad_covar)) {
            value = data$dyad_covar[[covar_lab]]
        } else {
            value = data$actor_covar[[covar_lab]]
        }

        if (covariate$is_within) {
            t = covariate$t_or_h
            D_within[t, var_name] = D_within[t, var_name] + 1
            if (sum(dim(value)) == 2*n) { #dyadic covar density, reciprocity,
                dim(value) <- c(n, n, 1)
                dimnames(value)[[3]] <- covar_name
                covariates_stan[[var_name]] <- abind::abind(covariates_stan[[var_name]], value, along = 3)
                #update prior

            } else { # actor covar  sender, receiver,
                new_names <- c(colnames(covariates_stan[[var_name]]), covar_name)
                covariates_stan[[var_name]] <- cbind(covariates_stan[[var_name]], value)
                colnames(covariates_stan[[var_name]]) <- new_names
            }
            
        } else { #cross_density, cross_reciprocity
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
                            cross_rho_covariates=covariates_stan$cross_reciprocity
                            )
    #print(dim(covariates_data$mu_covariates))
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
    obs_idx = which(!is.na(y)) # indices of observed dyads
    y_obs = y[obs_idx] # observed dyads

    network_data = list(n=n, N=N, T=t, H=H, layer_pairs=pairs,y_obs=y_obs, K=K, N_obs = length(y_obs),obs_idx = obs_idx)

    return(network_data)
}

#IP
format_prior <- function(baseline_prior, covariates, data) {
    dict_var_name_to_prior_name <- list("density" = c("mu_covariates_mean_prior", "mu_covariates_sd_prior"), 
                                        "reciprocity" = c("rho_covariates_mean_prior", "rho_covariates_sd_prior"),
                                        "sender" = c("alpha_covariates_mean_prior", "alpha_covariates_sd_prior"),
                                        "receiver" = c("beta_covariates_mean_prior", "beta_covariates_sd_prior"),
                                        "cross_density" = c("cross_mu_covariates_mean_prior", "cross_mu_covariates_sd_prior"),
                                        "cross_reciprocity" = c("cross_rho_covariates_mean_prior", "cross_rho_covariates_sd_prior"))
    prior <- list(mu_covariates_mean_prior = array(0, 0),
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

    for (covar_name in names(covariates)) {
        covariate = covariates[[covar_name]]
        var_name = covariate$var_name # sender, receiver, density, reciprocity, cross_density, cross_reciprocity
        covar_lab = covariate$covar_lab
        if (covar_lab %in% names(data$dyad_covar)) {
            value = data$dyad_covar[[covar_lab]]
        } else {
            value = data$actor_covar[[covar_lab]]
        }
        
        #update prior
        prior_mean_name = dict_var_name_to_prior_name[[var_name]][1]
        prior_sd_name = dict_var_name_to_prior_name[[var_name]][2]

        prior_mean_val = ifelse(is.null(covariate$mean_prior), 0, covariate$mean_prior)
        prior_sd_val = ifelse(is.null(covariate$sd_prior), 10/sd(value), covariate$sd_prior)

        new_names <- c(names(prior[[prior_mean_name]]), covar_name)
        prior[[prior_mean_name]] <- array(append(prior[[prior_mean_name]], prior_mean_val))
        prior[[prior_sd_name]] <- array(append(prior[[prior_sd_name]], prior_sd_val))
        #print(dim(prior[[prior_mean_name]]))
        names(prior[[prior_mean_name]]) <- new_names
        names(prior[[prior_sd_name]]) <- new_names

    }
    return(append(baseline_prior, prior))
}

#' The main function to convert user network input to Stan data
#' 
#' @param M list of the observed networks in adjacency matrix form
#' @param covariates intermediate form of the covariates: covariate <- list(is_within = {T, F}, t_or_h, var_name, value)
#' @return a list in the format for Stan code

create_stan_data <- function(model_obj) {
    stan_data_net = format_network(model_obj$data$dep_net)
    stan_data_covar = format_covariate(model_obj$t, model_obj$H, model_obj$n, model_obj$model$covar, data = model_obj$data)
    stan_data_prior = format_prior(model_obj$model$prior, model_obj$model$covar, model_obj$data)
    stan_data = c(stan_data_net, stan_data_covar, stan_data_prior)
    model_obj$model$full_prior = stan_data_prior
    model_obj$fit_res$stan_data = stan_data
    return(model_obj)
}


default_prior_empty_mdl <- function(dep_lab, pair_lab) {
    t = length(dep_lab)
    H = t*(t - 1)/2
    stopifnot(H == length(pair_lab))
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


    return(prior)
}

