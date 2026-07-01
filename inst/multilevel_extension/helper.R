flatten_dyad_covar <- function(covars) {
    
    flatten_dyad_covar_single <- function(X){
        res <- c()
        n = nrow(X)
        for (i in 1:n) {
            for (j in 1:n) {
                if (i != j) {
                res <- c(res,X[i,j])
                }
            }
        }
        return(res)
    }

    p = dim(covars)[3]
    n = dim(covars)[1]
    res <- lapply(1:p, function(i) flatten_dyad_covar_single(covars[,,i]))
    res <- do.call(rbind, res)
    return(res)

}

generate_stan_data_list <- function(L, n_l_seq, L_networks, layer_covar = FALSE) {
   stan_data_list <- lapply(1:L, function(l) {
    n_l = n_l_seq[l]
    dyad_covar <- list("dyad_covar1" = matrix(runif(n_l*n_l), nrow=n_l, ncol=n_l), 
                       "dyad_covar2" = matrix(runif(n_l*n_l), nrow=n_l, ncol=n_l))
    actor_covar <- matrix(runif(n_l*2), nrow=n_l, ncol=2)
    colnames(actor_covar) <- c("actor_covar1", "actor_covar2")
    actor_covar <- as.data.frame(actor_covar)

    dep_net <- L_networks[[l]]
    m_empty <- Mp2Model(dep_net, dyad_covar = dyad_covar, actor_covar = actor_covar)
    if (layer_covar) {
        m_empty <- update_covar(m_empty, layer_1 = "network1", density = c("dyad_covar1", "dyad_covar2"), sender = "actor_covar1", receiver = "actor_covar2")
        m_empty <- update_covar(m_empty, layer_1 = "network2", density = c("dyad_covar1", "dyad_covar2"), sender = "actor_covar1", receiver = "actor_covar2")
        m_empty <- update_covar(m_empty, layer_1= "network1", layer_2 = "network2", cross_density = c("dyad_covar1", "dyad_covar2"), cross_reciprocity = "dyad_covar1")
    }
    create_stan_data(m_empty)$fit_res$stan_data
    })
    return(stan_data_list)
}

create_stan_data_per_class <- function(class, data_path, layer_1 = "sts", layer_2 = "tst", layer_covar = FALSE, WAVE=4, covar_df = summary_df) {
    load(file.path(data_path, paste0(class,"_bully1_wave_",WAVE,"_cleaned.RData")))
    actor_covar_name = c("is_female", "ses")
    dep_net = list()
    dep_net[[layer_1]] = cleaned_data[[layer_1]]
    dep_net[[layer_2]] = cleaned_data[[layer_2]]
    actor_covar <- summary_df %>% filter(class_id == class, wave.x == WAVE) %>%
        select(all_of(actor_covar_name))
    dyad_covar <- list()
    # dyad_covar[["same_gender"]] <- outer(actor_covar[["is_female"]], actor_covar[["is_female"]], FUN = function(x,y) as.numeric(x==y))

 
    m_2 <- Mp2Model(dep_net, dyad_covar = dyad_covar, actor_covar = actor_covar)
    if (layer_covar) {
        m_2 <- update_covar(m_2, layer_1 = layer_1, receiver = actor_covar_name, sender = actor_covar_name)
        m_2 <- update_covar(m_2, layer_1 = layer_2, receiver = actor_covar_name, sender = actor_covar_name)
        # m_2 <- update_covar(m_2, layer_1 = layer_1, layer_2 = layer_2, cross_reciprocity = "same_gender")
        # m_2 <- update_covar(m_2, layer_1 = layer_1,  density = "same_gender")
        # m_2 <- update_covar(m_2, layer_1 = layer_2,  density = "same_gender")
    }
    stan_data <- create_stan_data(m_2)$fit_res$stan_data
    return(stan_data)
}


create_classlevel_covar <- function(class_id, var_name, summary_df) {
    covar <- mean(summary_df[summary_df$class_id == as.numeric(class_id),][[var_name]], na.rm = TRUE)
    return(covar)
}

generate_multilevel_stan_data <- function(stan_data_list){
    stan_data <- stan_data_list[[1]]
    stan_data$L = length(stan_data_list)
    stan_data$n <- sapply(stan_data_list, function(x) x$n)
    stan_data$N <- sapply(stan_data_list, function(x) x$N)
    stan_data$y_obs <- unlist(sapply(stan_data_list, function(x) x$y_obs))
    stan_data$obs_idx <- 1:length(stan_data$y_obs) #not used currently, assuming all observations are observed
    stan_data$network_sim = FALSE
    stan_data$prior_sim = FALSE


    return(stan_data)
}

correct_dim <- function(x) {
    if (length(x) != 0) {
        dim(x) <- c(length(x))
    }
    return(x)
}

add_group_covar <- function(stan_data, group_covariates=NULL, D_group_within=NULL, D_group_cross=NULL, mu_group_covariates_idx=numeric(0), rho_group_covariates_idx=numeric(0), cross_mu_group_covariates_idx=numeric(0), cross_rho_group_covariates_idx=numeric(0)) {
    create_covar_sd_prior <- function(covar_idx) {
        if (length(covar_idx) > 0) {
            res <- sapply(covar_idx, function(x) 10/sd(stan_data$group_covariates[,x]))
            res <- correct_dim(res)
            return(res)
        } else {
            return(numeric(0))
        }
    }


    t = stan_data$T
    H = stan_data$H
    L = stan_data$L
    if (is.null(group_covariates)) {
        group_covariates <- matrix(, nrow=L, ncol=0)
    }
    stan_data$group_covariates <- group_covariates
    stan_data$p_group <- ncol(stan_data$group_covariates)

    if (is.null(D_group_within)) {
        D_group_within <- matrix(0, nrow=2, ncol=t)
    }
    stan_data$D_group_within <- D_group_within

    if (is.null(D_group_cross)) {
        D_group_cross <- matrix(0, nrow=2, ncol=H)
    }
    stan_data$D_group_cross <- D_group_cross

    stan_data$mu_group_covariates_idx <- correct_dim(mu_group_covariates_idx)
    stan_data$rho_group_covariates_idx <- correct_dim(rho_group_covariates_idx)
    stan_data$cross_mu_group_covariates_idx <- correct_dim(cross_mu_group_covariates_idx)
    stan_data$cross_rho_group_covariates_idx <- correct_dim(cross_rho_group_covariates_idx)

    stan_data$mu_group_covariates_mean_prior <- correct_dim(rep(0, sum(stan_data$D_group_within[1,])))
    stan_data$rho_group_covariates_mean_prior <- correct_dim(rep(0, sum(stan_data$D_group_within[2,])))
    stan_data$cross_mu_group_covariates_mean_prior <- correct_dim(rep(0, sum(stan_data$D_group_cross[1,])))
    stan_data$cross_rho_group_covariates_mean_prior <- correct_dim(rep(0, sum(stan_data$D_group_cross[2,])))

    stan_data$mu_group_covariates_sd_prior <- create_covar_sd_prior(stan_data$mu_group_covariates_idx)
    stan_data$rho_group_covariates_sd_prior <- create_covar_sd_prior(stan_data$rho_group_covariates_idx)
    stan_data$cross_mu_group_covariates_sd_prior <- create_covar_sd_prior(stan_data$cross_mu_group_covariates_idx)
    stan_data$cross_rho_group_covariates_sd_prior <- create_covar_sd_prior(stan_data$cross_rho_group_covariates_idx)
    return(stan_data)
}

add_layer_covar <- function(stan_data, stan_data_list) {
    create_actor_covar <- function(total_num_covar, covar_name) {
        if (total_num_covar > 0) {
            res <- do.call(rbind, lapply(stan_data_list, function(x) x[[covar_name]]))
        } else {
            res <- matrix(, ncol = 0, nrow = sum(stan_data$n))
        }
        return(res)
    } 

    create_dyad_covar <- function(total_num_covar, covar_name) {
        if (total_num_covar > 0) {
            res <- lapply(stan_data_list, function(x) flatten_dyad_covar(x[[covar_name]]))
            res <- t(do.call(cbind, res))
        } else {
            res <- matrix(, ncol = 0, nrow = sum(stan_data$N_covar))
        }
        return(res)
    }

    aggregate_layer_covar_sd_prior <- function(covar_sd_name) {
        if (length(stan_data[[covar_sd_name]]) > 0) {
            res <- do.call(rbind, lapply(stan_data_list, function(x) x[[covar_sd_name]]))
            res <- correct_dim(colMeans(res))
        } else {
            res <- numeric(0)
        }

        return(res)

    }


    stan_data$alpha_covariates <- create_actor_covar(sum(stan_data$D_within[,3]), "alpha_covariates")
    stan_data$beta_covariates <- create_actor_covar(sum(stan_data$D_within[,4]), "beta_covariates")

    #combine all the dyad covariates:
    stan_data$N_covar = stan_data$n*(stan_data$n - 1) 
    stan_data$mu_covariates <- create_dyad_covar(sum(stan_data$D_within[,1]), "mu_covariates")
    stan_data$rho_covariates <- create_dyad_covar(sum(stan_data$D_within[,2]), "rho_covariates")
    stan_data$cross_mu_covariates <- create_dyad_covar(sum(stan_data$D_cross[,1]), "cross_mu_covariates")
    stan_data$cross_rho_covariates <- create_dyad_covar(sum(stan_data$D_cross[,2]), "cross_rho_covariates")

    #edit the layer covariate sd prior:
    stan_data$mu_covariates_sd_prior <- aggregate_layer_covar_sd_prior("mu_covariates_sd_prior")
    stan_data$rho_covariates_sd_prior <- aggregate_layer_covar_sd_prior("rho_covariates_sd_prior")
    stan_data$cross_mu_covariates_sd_prior <- aggregate_layer_covar_sd_prior("cross_mu_covariates_sd_prior")
    stan_data$cross_rho_covariates_sd_prior <- aggregate_layer_covar_sd_prior("cross_rho_covariates_sd_prior")
    stan_data$alpha_covariates_sd_prior <- aggregate_layer_covar_sd_prior("alpha_covariates_sd_prior")
    stan_data$beta_covariates_sd_prior <- aggregate_layer_covar_sd_prior("beta_covariates_sd_prior")

    return(stan_data)
}

# test_that(" find_dyadic_covar_l_idx works", {
    ##testing the find_dyadic_covar_l_idx function 
    # rstan::expose_stan_functions("/home/annihong/projects/multip2/inst/stan/dev_model.stan")

#    res <- c()
#     for (l in 1:L) {
#         for (i in 1:n_l_seq[l]) {
#             for (j in 1:n_l_seq[l]) {
#                 if (i != j) {
#                 res <- c(res, L_networks[[l]][[1]][i,j])
#                 }
#             }
#         }
#     }

#     for (l in 1:L) {
#         for (i in 1:n_l_seq[l]) {
#             for (j in 1:n_l_seq[l]) {
#                 if (i != j) {
#                 group_l_idx <- find_start_end(stan_data$N_covar, l)
#                 dyad_l_idx <- find_dyadic_covar_l_idx(n_l_seq[l], i, j)
#                 dyad_idx <- group_l_idx[1] + dyad_l_idx - 1
#                     if (res[dyad_idx] != L_networks[[l]][[1]][i,j]) {
#                         stop("Error in find_dyadic_covar_l_idx function")
#                     }
#                 }
#             }
#         }
#     }
# })