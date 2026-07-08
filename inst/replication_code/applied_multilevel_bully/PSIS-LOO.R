###############################################################################
# PSIS-LOO comparison of the empty vs group_only multilevel multiplex p2
# models (bully data, companion to model_fitting.R)
#
# No model is refit: rstan::gqs() re-runs ONLY the generated quantities block
# (with loo_est = 1) over the posterior draws saved by model_fitting.R.
# gqs() recomputes transformed parameters (x_beta) from the saved parameter
# draws, so it does not matter that x_beta was excluded when the fits were
# saved.
#
# NOTE: bully1_2000_full.Rds (2025-09-16) predates the 2026-02-16 regeneration
# of the cleaned data and was fit on buggy networks for classes 14100, 14200,
# 15100, 15200 -- it must be refit before it can join this comparison.
###############################################################################

devtools::load_all()
library(rstan)
library(dplyr)
library(loo)
options(mc.cores = 4)
rstan::rstan_options(auto_write = TRUE)
source("./inst/multilevel_extension/helper.R")

data_path <- "/home/annihong/projects/multiplex-social-networks/stan/coding/gossip_perception/Data/Raw/GossipBully"
fit_dir <- "/home/annihong/projects/multiplex-social-networks/multilevelp2/code"

###### bully data (identical to model_fitting.R) ##########

filtered_df <- read.csv("data-raw/bully/filtered_df_no_missing_wave_4_bully1.csv")
class_ids <- filtered_df$class_id
summary_df <- read.csv("/home/annihong/projects/multip2/data-raw/gossip/summary.csv")
summary_df <- summary_df %>% rename(class_id = class) %>% filter(class_id %in% class_ids, wave == 4)
summary_df <- summary_df %>% left_join(filtered_df, by = c("class_id")) %>% mutate(is_female = ifelse(gender == 2, 1, 0))

L <- length(class_ids)
t <- 2
H <- 1

# Rebuild the stan_data exactly as model_fitting.R does (R-side data prep only,
# no sampling).
build_stan_data <- function(layer_covar, group_covar) {
    stan_data_list <- lapply(class_ids, create_stan_data_per_class,
                             layer_covar = layer_covar, data_path = data_path,
                             layer_1 = "bully1", layer_2 = "bully1_rev", covar_df = summary_df)
    stan_data <- generate_multilevel_stan_data(stan_data_list)
    stan_data <- add_layer_covar(stan_data, stan_data_list)
    if (group_covar) {
        ses_percent <- sapply(class_ids, create_classlevel_covar, var_name = "ses", summary_df = summary_df)
        class_size <- sapply(class_ids, create_classlevel_covar, var_name = "num_present", summary_df = summary_df)
        group_covariates <- matrix(c(class_size, ses_percent), nrow = L, ncol = 2)
        D_group_within <- matrix(c(2, 0, 2, 0), nrow = 2, ncol = t)
        D_group_cross <- matrix(c(0, 0), nrow = 2, ncol = H)
        stan_data <- add_group_covar(stan_data, group_covariates, D_group_within, D_group_cross,
                                     mu_group_covariates_idx = c(1, 2, 1, 2))
    } else {
        stan_data <- add_group_covar(stan_data)
    }
    return(stan_data)
}

# Flags so the generated quantities block computes log_lik and skips the
# posterior network draws (y_tilde).
set_gq_flags <- function(stan_data) {
    stan_data$prior_sim <- 0
    stan_data$network_sim <- 0
    stan_data$loo_est <- 1
    return(stan_data)
}

# The saved fits must come from the same model version / covariate setup as the
# stan_data we feed to gqs(): the random-effect dimensions have to line up.
check_fit_data_compatibility <- function(fit, stan_data, label) {
    expected_z_within <- c(stan_data$T, 2 + sum(stan_data$D_within[1, ]), stan_data$L)
    actual_z_within <- fit@par_dims$z_within
    if (!identical(as.numeric(actual_z_within), as.numeric(expected_z_within))) {
        stop(label, ": z_within dims in saved fit (", paste(actual_z_within, collapse = "x"),
             ") do not match stan_data (", paste(expected_z_within, collapse = "x"),
             "); the fit comes from a different model version or covariate setup.")
    }
    expected_group_coef <- sum(stan_data$D_group_within[1, ])
    actual_group_coef <- fit@par_dims$mu_group_coef
    if (!identical(as.numeric(actual_group_coef), as.numeric(expected_group_coef))) {
        stop(label, ": mu_group_coef dims in saved fit do not match the group covariate setup.")
    }
    if (sum(stan_data$N) != length(stan_data$y_obs)) {
        stop(label, ": sum(N) does not match length(y_obs).")
    }
}

# Re-run only the generated quantities over the saved posterior draws and
# return the PSIS-LOO estimate.
run_psis_loo <- function(fit_file, stan_data, gq_model, label) {
    message(">>> ", label, ": loading fit ", fit_file)
    fit <- readRDS(fit_file)
    check_fit_data_compatibility(fit, stan_data, label)

    n_chains <- fit@sim$chains
    n_draws_per_chain <- (fit@sim$iter - fit@sim$warmup) / fit@sim$thin
    draws <- as.matrix(fit) # stacks chains: all draws of chain 1, then chain 2, ...
    rm(fit); gc()

    message(">>> ", label, ": running standalone generated quantities on ",
            nrow(draws), " draws")
    gq <- rstan::gqs(gq_model, data = stan_data, draws = draws)
    log_lik <- loo::extract_log_lik(gq, merge_chains = TRUE) # draws x sum(N), same row order as `draws`
    rm(gq, draws); gc()

    chain_id <- rep(1:n_chains, each = n_draws_per_chain)
    r_eff <- loo::relative_eff(exp(log_lik), chain_id = chain_id, cores = 4)
    loo_res <- loo::loo(log_lik, r_eff = r_eff, cores = 4)
    return(loo_res)
}

###### stan data ##########

stan_data_empty <- set_gq_flags(build_stan_data(layer_covar = FALSE, group_covar = FALSE))
stan_data_group_only <- set_gq_flags(build_stan_data(layer_covar = FALSE, group_covar = TRUE))
stan_data_full <- set_gq_flags(build_stan_data(layer_covar = TRUE, group_covar = TRUE))

# both models must be evaluated on the same outcomes for LOO to be comparable
stopifnot(identical(stan_data_empty$y_obs, stan_data_group_only$y_obs))
stopifnot(identical(stan_data_empty$y_obs, stan_data_full$y_obs))


###### PSIS-LOO ##########

# a compile (if not cached by auto_write), not a fit
gq_model <- rstan::stan_model("./inst/stan/multilevel_multiplex_p2_uncentered.stan")

# loo_empty <- run_psis_loo(file.path(fit_dir, "bully1_2000_empty.Rds"), stan_data_empty, gq_model, "empty")
# loo_group_only <- run_psis_loo(file.path(fit_dir, "bully1_2000group_only.Rds"), stan_data_group_only, gq_model, "group_only")


loo_full <- run_psis_loo(file.path(fit_dir, "bully1_2000full.Rds"), stan_data_full, gq_model, "full")

# saveRDS(loo_empty, file.path(fit_dir, "loo_bully1_2000_empty.Rds"))
# saveRDS(loo_full, file.path(fit_dir, "loo_bully1_2000_full.Rds"))
readRDS(file.path(fit_dir, "loo_bully1_2000_empty.Rds")) -> loo_empty
readRDS(file.path(fit_dir, "loo_bully1_2000_full.Rds")) -> loo_full
readRDS(file.path(fit_dir, "loo_bully1_2000_group_only.Rds")) -> loo_group_only


###### comparison ##########
print("empty model:")
print(loo_empty)
print("group-only model:")
print(loo_group_only)
print("full model:")
print(loo_full)

comparison <- loo::loo_compare(list(empty = loo_empty, group_only = loo_group_only, full = loo_full))
print(comparison, simplify = FALSE)
print(comparison, simplify = TRUE)
