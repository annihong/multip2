#load("../../results/sim_res_hydra6_num_sim_9_out_of_50.RData")

library(multip2)
path <- "/Users/annihong/Documents/Rprojects/simres/"
simulations <- list.files(paste0(path, "analysis_res/"))
simulations <- simulations[grep("test_", simulations)]
length(simulations)
#simulations <- grep("1000.Rds$", simulations, value = TRUE)
MODE = "PS"
#simulations <- simulations[grep("inv_gamma_hydra_5_.*_out_of_150", simulations)]
#simulations <- simulations[grep("inv_gamma", simulations)]
library(ggplot2)
library(dplyr)
# PS_mu_prior <- read.csv(file = "/home/annihong/projects/simres/analysis_res/PS_mu_prior.csv")[,c(2,3)]
# v_PS_mu <- (sd(unlist(PS_mu_prior)))^2

sampled_params_list <- list()
posterior_draws_list <- list()




zscore_helper <- function(true_param, vector) {
    post_mean = mean(vector)
    post_sd = sd(vector)
    zscore = (post_mean-true_param)/post_sd
    #print(paste0("post_mean: ", post_mean, "post_sd", post_sd))
    return(zscore)
}

post_contract_helper <- function(to_rank, vector, variance = 100) {
    post_sd = sd(vector)
    post_contract = 1 - (post_sd^2)/variance
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
    post_draws <- as.matrix(draws)
    stats <- c()

    for (t in 1:length(prior_params)) {
        prior_param = prior_params[[t]]
        if (MODE == "PS" & identical(stat_func, post_contract_helper) & param == "mu") {
            n = nrow(sampled_params$C)
            A_idx = 1 + 2 * (t - 1)
            B_idx = 2 + 2 * (t - 1)
            v_mu = 100 
            v_A = sampled_params$Sigma[A_idx, A_idx]
            v_B = sampled_params$Sigma[B_idx, B_idx]
            cov_AB = sampled_params$Sigma[A_idx, B_idx]
            v_PS_mu = v_mu + (1/n) * (v_A + v_B + 2 * cov_AB)
            stat <- stat_func(prior_param, post_draws[,t], variance = v_PS_mu)
        } else {
            stat <- stat_func(prior_param, post_draws[,t])
        }
        stats <- c(stats, stat)
    }
    names(stats) <- paste(param, 1:length(prior_params), sep = "_")
    return(stats)
}

calc_stats_file <- function(simulation_file, solution){
    file_path <- paste0(path, "analysis_res/", simulation_file)
    #print(simulation_file)
    # e1 <- new.env(parent = baseenv())
    # load(file_path, envir = e1)

    sim_result <- readRDS(file_path)
    stan_fit <- sim_result$Mp2_fit$fit_res$stan_fit
    sampled_params <- sim_result$sampled_params
    if (MODE == "OG") {
        posterior_draws <- rstan::extract(stan_fit, c("mu", "rho", "cross_mu", "cross_rho"))
    } else if (MODE == "PS"){
        posterior_draws <- rstan::extract(stan_fit, c("mu", "rho", "cross_mu", "cross_rho", "A_bar", "B_bar"))
        PS_mu = posterior_draws$mu + posterior_draws$A_bar + posterior_draws$B_bar
        posterior_draws$mu <- PS_mu
        posterior_draws <- posterior_draws[c("mu", "rho", "cross_mu", "cross_rho")]
        t_total = ncol(PS_mu)
        A_bar_prior <- colMeans(sampled_params$C[,1 + 2*(1:t_total - 1)])
        B_bar_prior <- colMeans(sampled_params$C[,2 + 2*(1:t_total - 1)])
        sampled_params$mu <- sampled_params$mu + A_bar_prior + B_bar_prior 
    
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
    res$sampled_params <- sampled_params
    res$posterior_draws <- posterior_draws
    return(res)
}


res <- lapply(simulations, calc_stats_file)
rho_s <- sapply(1:length(res), function(i) res[[i]]$sampled_params$rho)
mu_s <- sapply(1:length(res), function(i) res[[i]]$sampled_params$mu)

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
save(res, file=paste0(path, "analysis_res/sim_analysis_df_full_", MODE, "_1.RData"))
load(paste0(path, "analysis_res/sim_analysis_df_full_", MODE, "_1.RData"))
saveRDS(list(params=sampled_params_list, draws=posterior_draws_list), file=paste0(path, "analysis_res/prior_n_posterior_draws.rds"))

rank_df <- res$rank_df
zscore_df=res$zscore_df
post_contract_df=res$post_contract_df

hist <- ggplot(tidyr::gather(rank_df)) + geom_histogram(aes(x=value/101, color = key), size = 1, breaks=seq(0,1,0.05)) +    
        facet_wrap(~key, ncol=6) + ylab("") + xlab("") +
        #scale_y_continuous(labels = scales::percent) + 
        theme_bw() +
        scale_color_manual(values=c(rep("#fb8500", 1),  rep("#ffb703", 1), rep("#219ebc", 2), rep("#5AB1BB", 2))) + theme(legend.position="none", axis.text.x = element_text(size = 5), axis.text.y = element_text(size = 5)) #

    ggsave(paste0(path, "plots/hist_corrected_", MODE, "_1.png"), plot=hist, width = 8, height = 2)

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

ggsave(paste0(path, "plots/model_sensitivity_", MODE, "_1.png"), plot=scatter, width = 8, height = 3)



# library(ggplot2)
# set.seed(1)
# x<-runif(100,-1,1)
# dd<-data.frame(x)

# ks.test(rank_df$rho_2,"punif",1,101)

# rank_df <- res$rank_df
# rank_df <- rank_df[300:1000,]
# ed <- ecdf(rank_df$rho_2)
# maxdiffidx <- which.max(abs(ed(rank_df$rho_2)-punif(rank_df$rho_2,1,101)))
# maxdiffat <- rank_df$rho_2[maxdiffidx]

# p<-ggplot(aes(rho_2),data=rank_df)+stat_ecdf()+theme_bw()+stat_function(fun=punif,args=list(1,101))
# p<-p+labs(title="ECDF and theoretical CDF")+geom_vline(xintercept=maxdiffat, lty=2)
# p

