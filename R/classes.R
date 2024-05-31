

#' Format covariates data for stan estimation
#' 
#' @param outcome A character vector of the names of the layers of networks in the network_data as the outcome multiplex network
#' @param network_data list of p + t matrices of dim n x n : network outcomes and network covariates, p is the number of network covariates and t is the number of layers in the multiplex network outcome. 
#' @param actor_data n x p data frame: p covariates of length of the number of actors (n) 
#' @return a Mp2Model object 
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


#create a default empty model. 
Mp2Model <- function(dep_net, dyad_covar=NULL, actor_covar=NULL, ...) {
    # stopifnot(is_valid(dep_net, type = "dep_net"))
    # stopifnot(is_valid(dyad_covar, type = "dyad_covar"))
    # stopifnot(is_valid(actor_covar, type = "actor_covar"))

    t = length(dep_net)
    H = t*(t - 1)/2
    n = nrow(dep_net[[1]])
    dep_lab = names(dep_net)
    pair_lab = get_pair_names(dep_lab)

    covar = list() #start with an empty model
    prior <- default_prior_empty_mdl(t, H, dep_lab, pair_lab)

    newMp2Model = structure(list(
        t = t, H = H, n = n,
        dep_lab =dep_lab, pair_lab = pair_lab,
        dyad_covar_lab = names(dyad_covar), actor_covar_lab = names(actor_covar),
        data = list(dep_net = dep_net, dyad_covar = dyad_covar, actor_covar = actor_covar),
        model = list(covar=covar, prior=prior),
        fit_res = list(stan_data = NULL, stan_fit = NULL, summary = NULL, par_labels = NULL)),
        class = "Mp2Model")

    return(newMp2Model)
}

#' Add covariates to the model
add_covar <- function(model_obj,...) {
  UseMethod("add_covar")
}

# param = c("density", "reciprocity", "sender", "receiver"), c("cross_density", "cross_reciprocity")
#list(var_name=var_name, t_or_h=t_or_h, is_within=is_within, value=value)

add_covar.Mp2Model <- function(model_obj, covar_lab, param, layer_lab=NULL, mean_prior=NULL, sd_prior=NULL ...) {
    create_covar <- function(var_name, t_or_h, is_within, value) {
        return(list(var_name=var_name, t_or_h=t_or_h, is_within=is_within, value=value))
    }
    within_param <- c("density", "reciprocity", "sender", "receiver")
    cross_param <- c("cross_density", "cross_reciprocity")
    dyad_param <- c("density", "reciprocity", "cross_density", "cross_reciprocity")
    actor_param <- c("sender", "receiver")

    if (param %in% dyad_param) {
        value = model_obj$data$dyad_covar[[covar_lab]]
    } else if (param %in% actor_param) {
        value = model_obj$data$actor_covar[[covar_lab]]
    } else {
        stop("Invalid param")
    }
    if (is.null(value)) {
        stop("Covariate not found")
    }

    #names(res) = paste(var_name, covar, outcome, sep="_")
    if (param %in% within_param) {
        is_within = TRUE
        if (is.null(layer_lab)) {
            idx_seq = 1:model_obj$t
            layer_lab = model_obj$dep_lab
        }else{
            idx_seq = which(model_obj$dep_lab == layer_lab)
        }
    } else if (param %in% cross_param) {
        is_within = FALSE
        if (is.null(layer_lab)) {
            idx_seq = 1:model_obj$H
            layer_lab = model_obj$pair_lab
        }else{
            idx_seq = which(model_obj$pair_lab == layer_lab)
        }
    } else {
        stop("Invalid parameter")
    }

    res = lapply(idx_seq, create_covar, var_name=param, is_within=is_within, value=value)
    names(res) = paste(param, covar_lab, layer_lab, sep="_")
    model_obj$model$covar = append(model_obj$model$covar, res)

    return(model_obj)
}


#' fit function
#' 
#' This function performs the estimation of the multiplex P2 model using Stan.
#' 
#' @param model_obj An object of class "Mp2Model"
#' @return An object of class "Mp2Model" with the fitted Stan model and parameter labels
#' @export
fit <- function(model_obj, ...) {
  UseMethod("fit")
}

#' fit function for "Mp2Model" class
#'
#' @param model_obj An object of class "Mp2Model"
#' @param ... Other arguments passed to or from other methods
#'
#' @export fit.Mp2Model
#' @export
fit.Mp2Model <- function(model_obj, chains = 4, iter = 200, warmup = floor(iter/2),
                     thin = 1, seed = sample.int(.Machine$integer.max, 1),
                     stan_file = "multiplex_p2_low_mem.stan", prior_sim = FALSE,mc.cores = parallel::detectCores(), auto_write = TRUE,...) {

    stan_data <- create_stan_data(model_obj)
    stan_data$prior_sim = prior_sim

    options(mc.cores =mc.cores)
    rstan::rstan_options(auto_write = auto_write)
    
    fpath <- system.file("stan", stan_file, package="multiplexP2")
    p2_fit <- rstan::stan(
                    file=fpath, 
                    data = stan_data,
                    chains = chains,
                    iter = iter,
                    warmup = warmup,
                    thin = thin,
                    seed = seed
                    )
    model_obj$fit_res$stan_data <- stan_data
    model_obj$fit_res$stan_fit <- p2_fit
    s <- create_summary(model_obj)
    model_obj$fit_res$summary <- s$summary
    model_obj$fit_res$par_labels <- s$par_labels
    return(model_obj)
}

#' Create summary for a model object
#'
#' This function creates a summary for a given model object. It extracts fixed and random effects summaries
#' and returns them along with parameter labels.
#'
#' @param model_obj The model object for which the summary is to be created.
#' @return A list containing the fixed and random effects summaries, along with parameter labels.
#' @export
create_summary <- function(model_obj) {
    pair_lab= model_obj$pair_lab
    dep_lab = model_obj$dep_lab
    fixed_summary <- make_fixed_summary(p2_fit, stan_data, dep_lab, pair_lab)
    random_summary <- make_random_summary(p2_fit, dep_lab)         
    summary <- list("fixed" = fixed_summary$summary, "random" = random_summary$summary)
    par_labels <- rbind(fixed_summary$par_labels, random_summary$par_labels)
    rownames(par_labels) <- NULL
    return(list(summary = summary, par_labels = par_labels))
}

#' get draws function
#' 
#' This function obtains the draws from the fitted Stan model of the Mp2Model class.
#' 
#' @param model_obj An object of class "Mp2Model"
#' @return draws in the speficied format
#' @export
extract_network_draws <- function(model_obj, ...) {
  UseMethod("extract_network_draws")
}

#' Extract simulated network outcome from the (prior) posterior of a fitted stan object of the Mp2Model class
#' 
#' @param model_obj  An object of class "Mp2Model"
#' @param sim_num integer: number of simulations to extract, counting from the tail of the posterior draws
#' @param network_type Type of network to return. Options are "adj" for adjacency matrix (default), "igraph" for igraph object, or "network" for network object, "dyad" for dyad form.
#' @return A list of length `sim_num`. Each element of the list is another list of `t` networks. These networks represent the (prior) posterior draws of the simulated network outcome, in the specified format.
#' @export
extract_network_draws <- function(model_obj, sim_num, network_type = "adj") {
    n <- attr(model_obj, "stan_data")$n
    t <- attr(model_obj, "stan_data")$T
    fit <- model_obj$stan_fit
    network_draws <- extract_draws(fit, "y_tilde")
    res <- tail.matrix(network_draws, sim_num)
    if (network_type != "dyad") {
        res <- dyads_to_matrix_list(dyad_df = res, n = n, t = t, network_type = network_type)
    }
    return(res)
}


#' @export
summary.Mp2Model <- function(model_obj, ...) {
    s <- model_obj$summary
    print(s)
    
}


#' @export
print.Mp2Model <- function(model_obj, ...) {
    cat(attr(model_obj, "outcome"))
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
#' @param stan_data list: the "stan_data" attribute of the Mp2Model object
#' @param outcome vector: the "outcome" attribute of the Mp2Model object, names of the layers of the multiplex network
#' @param pairs vector: the "pair_names" attribute of the Mp2Model object, names of the pairs of layers of the multiplex network
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
#' @param outcome vector: the "outcome" attribute of the Mp2Model object, names of the layers of the multiplex network
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


