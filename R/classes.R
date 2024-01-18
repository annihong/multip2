
#' Format covariates data for stan estimation
#' 
#' @param outcome A character vector of the names of the layers of networks in the network_data as the outcome multiplex network
#' @param network_data list of p + t matrices of dim n x n : network outcomes and network covariates, p is the number of network covariates and t is the number of layers in the multiplex network outcome. 
#' @param actor_data n x p data frame: p covariates of length of the number of actors (n) 
#' @return a MultiP2Fit object 
#' @param senders_covar A character vector indicating the names of the covariates to include for the senders.
#' @param receivers_covar A character vector indicating the names of the covariates to include for the receivers.
#' @param density_covar A character vector indicating the names of the covariates to include for the density.
#' @param reciprocity_covar A character vector indicating the names of the covariates to include for the reciprocity.
#' @param cross_density_covar A character vector indicating the names of the covariates to include for the cross-layer density.
#' @param cross_reciprocity_covar A character vector indicating the names of the covariates to include for the cross-layer reciprocity.
#' @param custom_covars A list of custom covariates to include in the model.
#' @param chains An integer indicating the number of Markov chains to run.
#' @param iter An integer indicating the number of iterations per chain.
#' @param warmup An integer indicating the number of warmup iterations per chain.
#' @param thin An integer indicating the thinning rate for the chains.
#' @param seed An integer indicating the random seed for the chains.
#' @param stan_file A character string indicating the name of the stan file to use, default is "multiplex_p2.stan".
#' @return A list containing the fitted Stan model and parameter labels.
#' @export

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
                     thin = 1, seed = sample.int(.Machine$integer.max, 1),
                     stan_file = "multiplex_p2.stan", prior_sim = FALSE) {

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
        covariates = covariates_helper(t, H, covars, outcome, network_data, actor_data)
    } else {
        covariates = custom_covars
    }
    
    stan_covar = format_covariates(t, H, n, covariates)
    stan_data = append(stan_network_data, stan_covar)
    stan_data$prior_sim = prior_sim

    options(mc.cores = parallel::detectCores())
    rstan::rstan_options(auto_write = TRUE)
    
    fpath <- system.file("stan", stan_file, package="multiplexP2")
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
    pairs=unlist(get_pair_names(outcome)) 
    fixed_summary <- make_fixed_summary(p2_fit, stan_data, outcome, pairs)
    random_summary <- make_random_summary(p2_fit, outcome)         
    summary <- list("fixed" = fixed_summary$summary, "random" = random_summary$summary)
    par_labels <- rbind(fixed_summary$par_labels, random_summary$par_labels)
    rownames(par_labels) <- NULL
    newMultiP2Fit = structure(
                            list(stan_fit = p2_fit, summary = summary, par_labels = par_labels),
                            network_data=network_data,
                            actor_data=actor_data,
                            outcome=outcome,
                            pair_names=pairs,
                            covariates=covariates,
                            stan_data = stan_data,
                            class = "MultiP2Fit"
                            )

    return(newMultiP2Fit)
}


#' @export
summary.MultiP2Fit <- function(x, ...) {
    s <- x$summary
    print(s)
    
}


#' @export
print.MultiP2Fit <- function(x, ...) {
    cat(attr(x, "outcome"))
}

#' Extract parameter summary from the fitted stan object
#' 
#' @param fit rstan fit object: fitted stan object
#' @param pattern string: a regex pattern corresponding to the desired parameters, 
#' default is "^mu|^rho|^cross|fixed", which are all the fixed parameters
#' @return a matrix of the model output summary of the specified parameters
#' @export
extract_model_info <- function(fit, pattern = "^mu|^rho|^cross|fixed") {
  s <- rstan::summary(fit)$summary
  res <- s[grep(pattern,rownames(s)),1:10]
  return(res)
}

#' Extract and rename the summary of all the fixed parameters
#' 
#' @param fit rstan fit object: fitted stan object
#' @param stan_data list: the "stan_data" attribute of the MultiP2Fit object
#' @param outcome vector: the "outcome" attribute of the MultiP2Fit object, names of the layers of the multiplex network
#' @param pairs vector: the "pair_names" attribute of the MultiP2Fit object, names of the pairs of layers of the multiplex network
#' @return a matrix of the model output summary of the fixed parameters (baseline and covariates effects)
make_fixed_summary <- function(fit, stan_data, outcome, pairs) {

    rename_covar <- function(res, pattern, param_names, get_name_func) {
    network_covar_idx <- grep(pattern, rownames(res))
    network_covar_names <- lapply(param_names, get_name_func)
    names(network_covar_names) <- NULL
    newnames <- rownames(res)
    newnames[network_covar_idx] <- unlist(network_covar_names)
    rownames(res) <- newnames
    return(res)
    }

    pattern = "^mu|^rho|^cross|fixed"
    res = extract_model_info(fit=fit, pattern=pattern)
    old_names <- rownames(res)

    res <- rename_covar(res=res, pattern="^mu.*fixed|^rho.*fixed|^cross_mu.*fixed|^cross_rho.*fixed",
        param_names=c("mu_covariates", "rho_covariates", "cross_mu_covariates", "cross_rho_covariates"),
        get_name_func=function(covar) dimnames(stan_data[[covar]])[[3]]
                        )
    
    res <- rename_covar(res=res, pattern="^alpha.*fixed|^beta.*fixed",
    param_names=c("alpha_covariates", "beta_covariates"),
    get_name_func=function(covar) colnames(stan_data[[covar]])
                        )
    
    row_names <- rownames(res)
    # outcome <- attr(x, "outcome")
    # pairs <- attr(x, "pair_names")
    mu_idx <- grep("^mu\\[\\d\\]", row_names, perl = T)
    rho_idx <- grep("^rho\\[\\d\\]", row_names, perl = T)
    cross_mu_idx <- grep("^cross_mu\\[\\d\\]", row_names, perl = T)
    cross_rho_idx <- grep("^cross_rho\\[\\d\\]", row_names, perl = T)
    row_names[mu_idx] <- paste("density", outcome, sep = "_")
    row_names[rho_idx] <- paste("reciprocity", outcome, sep = "_")
    row_names[cross_mu_idx] <- paste("cross_density", pairs, sep = "_")
    row_names[cross_rho_idx] <- paste("cross_reciprocity", pairs, sep = "_")
    rownames(res) <- row_names
     P <- data.frame(
          Parameter = old_names,
          Label = row_names)
    return(list("summary" = res, "par_labels" = P))
}

#' Extract and rename the summary of all the random parameters
#' 
#' @param fit rstan fit object: fitted stan object
#' @param outcome vector: the "outcome" attribute of the MultiP2Fit object, names of the layers of the multiplex network
#' @return a matrix of the model output summary of the random parameters (the estimated variance-covariance matrix)
make_random_summary <- function(fit, outcome) {
    #extract all the entries of the varcov matrix from the p2 fit 
    sigma_raw <- extract_model_info(fit, pattern="Sigma")
    old_names <- rownames(sigma_raw)
    #names of the actor vars 
    new_names <- c(outer(c("sender:", "receiver:"),outcome, FUN=paste0))
    res <- sigma_raw
    for (name in rownames(sigma_raw)) {
        name_string = strsplit(name, "")[[1]]
        row = name_string[7]
        col = name_string[9]
        #remove all the repeated entries because the var-covar matrix is symmetrical 
        if (row > col) {
            res = res[-which(rownames(res) == name),]
        }
    }

    rename_sigma <- function(name) {
        name_string = strsplit(name, "")[[1]]
        row = name_string[7]
        col = name_string[9]
        new_name <- paste(new_names[as.numeric(row)], new_names[as.numeric(col)], sep = "_")
        return(new_name)
    }

    rownames(res) <- sapply(rownames(res), rename_sigma)
    P <- data.frame(
    Parameter = old_names,
    Label = sapply(old_names, rename_sigma))
    return(list("summary" = res, "par_labels" = P))
}

#' Extract draws from the (prior) posterior of a fitted stan object
#' 
#' @param fit rstan fit object: fitted stan object
#' @param parameter string: name of the desired parameters, y_tilde is the simulated outcome of the model
#' @return the (prior) posterior draws of the specified parameters
extract_draws <- function(fit, parameter) {
  s <- rstan::extract(fit)
  res <- s[[parameter]]
  return(res)
}

#' Extract simulated network outcome from the (prior) posterior of a fitted stan object
#' 
#' @param fit rstan fit object: fitted stan object
#' @param sim_num integer: number of simulations to extract, counting from the tail of the posterior draws
#' @param as_adjacency logical: whether to return the simulated network outcome as an adjacency matrix, default is TRUE
#' @return the (prior) posterior draws of the simulated network outcome
#' @export
extract_network_draws <- function(fit, sim_num, as_adjacency = TRUE) {
    network_draws <- extract_draws(fit, "y_tilde")
    res <- tail.matrix(network_draws, sim_num)
    return(res)
}
