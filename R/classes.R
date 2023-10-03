
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
        covariates = covariates_helper(t, H, covars, outcome, network_data, actor_data)
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
    pairs=unlist(get_pair_names(outcome))                
    summary <- list("fixed" = make_fixed_summary(p2_fit, stan_data, outcome, pairs), "random" = make_random_summary(p2_fit, outcome))
    newMultiP2Fit = structure(
                            list(stan_fit = p2_fit, summary = summary),
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

# networks = replicate(4, matrix(data=1:9, nrow=3), simplify = FALSE)
# names(networks) <-  c("network1", "network2", "network_covariate1",      "network_covariate2")
# actor_data = data.frame(actor_attr1=1:3, actor_attr2=1:3)


# m <- MultiP2Fit(outcome = c("network1", "network2"), 
# network_data = networks, 
# actor_data = actor_data, senders_covar=c("actor_attr1", "actor_attr2"), iter=100,
#             receivers_covar = c("actor_attr1", "actor_attr2"),
#             density_covar = c("network_covariate1","network_covariate2"),
#             reciprocity_covar = c("network_covariate1","network_covariate2"),
#             cross_density_covar = c("network_covariate1","network_covariate2"))


#' @export
summary.MultiP2Fit <- function(x) {
    s <- x$summary
    print(s)
    
}


#' @export
print.MultiP2Fit <- function(x) {
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
    return(res)
}

#' Extract and rename the summary of all the random parameters
#' 
#' @param fit rstan fit object: fitted stan object
#' @param outcome vector: the "outcome" attribute of the MultiP2Fit object, names of the layers of the multiplex network
#' @return a matrix of the model output summary of the random parameters (the estimated variance-covariance matrix)
make_random_summary <- function(fit, outcome) {
    #extract all the entries of the varcov matrix from the p2 fit 
    sigma_raw <- extract_model_info(fit, pattern="Sigma")
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
        # renaming the lower triangular entries but rename for easier interpretation
        } else {
            new_name <- paste(new_names[as.numeric(row)], new_names[as.numeric(col)], sep = "_")
            rownames(res)[rownames(res) == name] <- new_name
        }
    }
    return(res)
}

