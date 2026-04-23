# setwd("./multip2")
devtools::load_all()
library(rstan)
library(dplyr)
options(mc.cores = 4)
rstan::rstan_options(auto_write = TRUE)
source("./inst/multilevel_extension/helper.R")
data_path = "/home/annihong/projects/multiplex-social-networks/stan/coding/gossip_perception/Data/Raw/GossipBully"




######gossip data##########

filtered_df <- read.csv("data-raw/bully/filtered_df_no_missing_wave_4_bully1.csv")
class_ids = filtered_df$class_id
summary_df <- read.csv("/home/annihong/projects/multip2/data-raw/gossip/summary.csv")
summary_df <- summary_df %>% rename(class_id = class) %>% filter(class_id %in% class_ids, wave == 4) 
summary_df <- summary_df %>% left_join(filtered_df, by = c("class_id")) %>% mutate(is_female = ifelse(gender == 2, 1, 0))

# #######extension to L networks ######


L = length(class_ids)
t = 2
H = 1
layer_covar = FALSE
group_covar = TRUE



#########group covar###########

roma_percent <- sapply(class_ids, create_classlevel_covar, var_name = "roma", summary_df = summary_df)
ses_percent <- sapply(class_ids, create_classlevel_covar, var_name = "ses", summary_df = summary_df)
class_size <- sapply(class_ids, create_classlevel_covar, var_name = "num_present", summary_df = summary_df)

group_covariates <- matrix(c(class_size, ses_percent), nrow = L, ncol = 2)
D_group_within <- matrix(c(2,0,2,0), nrow = 2, ncol = t)
D_group_cross <- matrix(c(0,0), nrow = 2, ncol = H)

mu_group_covariates_idx <- c(1,2,1,2)
rho_group_covariates_idx <- numeric(0)
cross_mu_group_covariates_idx <- numeric(0)
cross_rho_group_covariates_idx <- numeric(0)


stan_data_list <- lapply(class_ids, create_stan_data_per_class, layer_covar = layer_covar, data_path = data_path, layer_1 = "bully1", layer_2 = "bully1_rev", covar_df = summary_df)


stan_data <- generate_multilevel_stan_data(stan_data_list)
stan_data <- add_layer_covar(stan_data, stan_data_list)
if (group_covar) {
    stan_data <- add_group_covar(stan_data, group_covariates, D_group_within, D_group_cross,
                              mu_group_covariates_idx, rho_group_covariates_idx,
                              cross_mu_group_covariates_idx, cross_rho_group_covariates_idx)
} else {
    stan_data <- add_group_covar(stan_data)
}

res <- rstan::stan(file = "./inst/stan/multilevel_multiplex_p2_uncentered.stan", data = stan_data, chains = 4, iter = 2000, par = c("x_beta"), include = FALSE)

if (layer_covar & group_covar) {
    model_type = "full"
} else if (!layer_covar & !group_covar) {
    model_type = "empty"
} else if(layer_covar & !group_covar) {
    model_type = "layer_only"
} else if(!layer_covar & group_covar) {
    model_type = "group_only"
}
saveRDS(res, file = glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/code/bully1_2000{model_type}.Rds"))
# saveRDS(stan_data, file = glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/code/stan_data_bully1_2000{ifelse(layer_covar, '_full', '_empty')}.Rds")


# # fit <- readRDS(file = glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/code/bully1_2000_empty.Rds"))
# fit <- readRDS(file = glue::glue("/home/annihong/projects/multiplex-social-networks/multilevelp2/code/bully1_2000_full.Rds"))
# fit@sim$iter          # total iterations per chain
# fit@sim$warmup        # warmup iterations                                                            
# fit@sim$iter - fit@sim$warmup  # sampling iterations  

