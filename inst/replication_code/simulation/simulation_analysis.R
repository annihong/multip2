#load("../../results/sim_res_hydra6_num_sim_9_out_of_50.RData")

library(multip2)
path <- "/home/annihong/projects/simres/"
simulations <- list.files(paste0(path, "sim_res/"))
length(simulations)
#simulations <- grep("1000.Rds$", simulations, value = TRUE)


#simulations <- simulations[grep("inv_gamma_hydra_5_.*_out_of_150", simulations)]
#simulations <- simulations[grep("inv_gamma", simulations)]
library(ggplot2)
library(dplyr)

zscore_helper <- function(true_param, vector) {
    post_mean = mean(vector)
    post_sd = sd(vector)
    zscore = (post_mean-true_param)/post_sd
    #print(paste0("post_mean: ", post_mean, "post_sd", post_sd))
    return(zscore)
}

post_contract_helper <- function(to_rank, vector) {
    post_sd = sd(vector)
    post_contract = 1 - (post_sd^2)/(10^2)
    return(post_contract)
}

rank_helper <- function(to_rank, vector){
    subsample <- sample(vector,100)
    # print(length(subsample))
    rank <- rank(c(to_rank, subsample))[1]
    return(rank)
}
# given the parameter name, 
calc_stats <- function(sampled_params, posterior_draws, param, stat_func){
    prior_params <- sampled_params[[param]]
    draws <- posterior_draws[[param]]
    # print(param)
    # print(nrow(draws))
    post_draws <- as.matrix(draws[sample(1:nrow(draws), 1000, replace=T),])
    stats <- c()
    for (i in 1:length(prior_params)) {
        prior_param = prior_params[[i]]
        stat <- stat_func(prior_param, post_draws[,i])
        stats <- c(stats, stat)
    }
    names(stats) <- paste(param, 1:length(prior_params), sep = "_")
    return(stats)
}

calc_stats_file <- function(simulation_file, solution){
    file_path <- paste0(path, "sim_res/", simulation_file)
    #print(simulation_file)
    # e1 <- new.env(parent = baseenv())
    # load(file_path, envir = e1)

    sim_result <- readRDS(file_path)
    stan_fit <- sim_result$Mp2_fit$fit_res$stan_fit
    sampled_params <- sim_result$sampled_params
    if (solution == "HC") {
        posterior_draws <- rstan::extract(stan_fit, c("mu", "rho", "cross_mu", "cross_rho"))
    } else if (solution == "PS"){
        posterior_draws <- rstan::extract(stan_fit, c("PS_mu", "rho", "cross_mu", "cross_rho"))
    }
    
    
    #posterior_draws <- sim_result$posterior_draws

    params <- names(sampled_params)[1:4] # no sigma
    names(posterior_draws) <- params # rename PS_mu as mu
    res <- list(ranks = c(), zscores = c(), posts = c())
    stat_funcs <- c(rank_helper, zscore_helper, post_contract_helper)
    for (param in params) {
        for (i in 1:length(res)) {
            res[[i]] = c(res[[i]], calc_stats(sampled_params, posterior_draws, param, stat_funcs[[i]]))
        }
    }
    return(res)
}

#  calc_z_score <- function(sampled_params, posterior_draws, param_name){
#         prior_param <- sampled_params[[param_name]]
#         true_param <- rep(0, length(prior_param))
#         post_draws <- posterior_draws[,
#         grep(paste0("^", param_name), colnames(posterior_draws))]
#         if (length(prior_param) > 1) {
#             ranks <- purrr::map2(prior_param, post_draws, rank_helper)
#         } else {
#         ranks <- rank_helper(prior_param, post_draws)
#         }
#         return(ranks)

#  }

res <- lapply(simulations, calc_stats_file, solution="PS")
rank_df <- data.frame()
zscore_df <- data.frame()
post_contract_df <- data.frame()
for (file_res in res) {
    rank_df <- rbind(rank_df,file_res$ranks)
    zscore_df <- rbind(zscore_df,file_res$zscores)
    post_contract_df <- rbind(post_contract_df, file_res$posts)
}

colnames(rank_df) <- colnames(zscore_df) <- colnames(post_contract_df) <- names(res[[1]]$ranks)
res=list(rank_df=rank_df, zscore_df=zscore_df, post_contract_df=post_contract_df)
save(res, file=paste0(path, "analysis_res/sim_analysis_df_full.RData"))
load(paste0(path, "analysis_res/sim_analysis_df_full.RData"))

# scale_fill_manual(values=c(rep("#fb8500", 2),  rep("#ffb703", 2), rep("#219ebc", 2), rep("#8ecae6", 2))) +

# for (col in colnames(rank_df)) {
#     rank_df[col] <- as.numeric(unlist(rank_df[col]))
# } , y = after_stat(count / sum(count))
#plot_rank <- 
hist <- ggplot(tidyr::gather(rank_df)) + geom_histogram(aes(x=value/101, color = key), size = 1, breaks=seq(0,1,0.1)) +    
        facet_wrap(~key, ncol=6) + ylab("") + xlab("") +
        #scale_y_continuous(labels = scales::percent) + 
        theme_bw() +
        scale_color_manual(values=c(rep("#fb8500", 1),  rep("#ffb703", 1), rep("#219ebc", 2), rep("#5AB1BB", 2))) + theme(legend.position="none", axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5)) #

    ggsave(paste0(path, "plots/hist_corrected.png"), plot=hist, width = 8, height = 2)

scatters <- list()
for (param_name in colnames(post_contract_df)) {
    df = data.frame(zscore=zscore_df[,param_name], post_contract=post_contract_df[,param_name])
    p <- ggplot(df) + geom_point(aes(x=post_contract, y=zscore)) + ylim(c(-5,5)) + ggtitle(param_name)
    scatters <- append(scatters, p)
    #ggsave(paste0(path, "plots/", param_name,".png"), plot=p)
}
zscore_long <- tidyr::gather(zscore_df, value = "zscore", key = "params")
scatter_df <- tidyr::gather(post_contract_df, value = "post_contract", key = "params")
scatter_df[,"z_score"] <- zscore_long$zscore

scatter <- ggplot(scatter_df,aes(x=post_contract,y=z_score))+ geom_point(aes(color=params), alpha = 0.4, size = 1) + facet_wrap( ~ params, ncol=6) + xlab("Posterior contraction") + ylab("Z-score") + 
theme_bw() + 
scale_color_manual(values=c(rep("#fb8500", 1),  rep("#ffb703", 1), rep("#219ebc", 2), rep("#5AB1BB", 2))) + theme(legend.position="none", axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5),
      axis.title = element_text(size = 15))

ggsave(paste0(path, "plots/model_sensitivity.png"), plot=scatter, width = 8, height = 3)

