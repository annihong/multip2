file_path <- "/home/annihong/projects/simres/fitted_models/"
seed = 102224

##LIBRARIES##
require(tidyverse)
library(multip2)
# require(ggmcmc)
# require(igraph)
# require(network)
# require(ergm)
# library(viridis)
# require(statnet.common)


##READ IN DATA##
results <- readRDS(paste0(file_path, "study_1_models.Rds"))



##LOAD IN THE STAN FIT OBJECTS##
fit1 <- results$m_1$fit_res$stan_fit
fit2 <- results$m_2$fit_res$stan_fit

##CONVERGENCE CHECK## NEEDS REVISION (figure out PS or HC, change the stan model used to estimation)

#S <- ggs(fit2, family = "PS|rho|Corr|fixed")
# S <- ggs(fit2, family = "PS")
# S <- ggs(fit2, family = "^Corr\[(\d+),(?=\d+\]$)(\d+)\]$")
# ggmcmc(S, file = paste0("diagnostics",".pdf"),  family = "PS|rho", param_page=6) #plot=c("traceplot", "running", "geweke"),
# ptrace <- ggs_traceplot(S, family = "^Corr") +
#       facet_wrap(~ Parameter, ncol = 6, scales = "free")
# p <- ggs_grb(S, family = "^mu") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
# p <- ggs_grb(S, family = "PS_mu") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
# ggs_geweke(S, family = "mu|rho") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
# ggs_running(S, family = "Sigma") + facet_wrap(~ Parameter, ncol = 6, scales = "free")
# ggs_density(ggs(radon$s.radon, par_labels=P, family="sigma"))
# ptrace <- ggs_traceplot(S, family = "mu")
# ggs_running(S, family = "^PS_mu")
# ggs_Rhat(S, family = "^mu")
# ggs_geweke(S, family = "^mu")
# ggsave(ptrace)

##MODEL ESTIMATES##

s1 <- summary.Mp2Model(results$m_1)
fit1 <- s1$fit_res$stan_fit

s2 <- summary.Mp2Model(results$m_2)
fit2 <- s2$fit_res$stan_fit


##### RANDOM EFFECTS ###### NEED TO FIX THIS!
t = results$m_1$t
multi_result_corr <- s2$correlation

# labs <- rownames(multi_result_corr)
# sender_labs <- grep("^sender", labs, value = T)
# receiver_labs <- grep("^receiver", labs, value = T)
labs <- c(paste("sender", c("PO", "SC", "PE"), sep =":"), paste("receiver", c("PO", "SC", "PE"), sep =":"))
corr_M <- matrix(, ncol = 2*t, nrow = 2*t, dimnames = list(labs, labs))

fill_corr_M <- function(corr_M, res) {
      res_M <- corr_M 
      for (col in colnames(corr_M)){
            for (row in rownames(corr_M)){
                  lab <- paste0(col, "_", row)
                  if (lab %in% names(res)) {
                        res_M[col, row] = res[lab]
                  } else {
                        lab <- paste0(row, "_", col)
                        res_M[col, row] = res[lab]
                  }
            }
      }
      return(res_M)
}


#sig_dat <- fill_corr_M(corr_M, multi_result_corr[, '2.5%'] * multi_result_corr[, '97.5%'] > 0 )
# val_dat <- fill_corr_M(corr_M, multi_result_corr[, "mean"])
# upper_dat <- fill_corr_M(corr_M, multi_result_corr[, '97.5%'])
# lower_dat <-fill_corr_M(corr_M, multi_result_corr[, '2.5%'])
df_corr <- expand.grid(Var1 = rownames(corr_M), Var2 = colnames(corr_M))
df_corr$value <- as.vector(fill_corr_M(corr_M, multi_result_corr[, "mean"]))
df_corr$upper <- as.vector(fill_corr_M(corr_M, multi_result_corr[, '97.5%']))
df_corr$lower <- as.vector(fill_corr_M(corr_M, multi_result_corr[, '2.5%']))
df_corr$CI_low <- paste0("(", round(df_corr$lower , 1))
df_corr$CI_high <- paste0(round(df_corr$upper, 1), ")")
df_corr$CI <- paste0(df_corr$CI_low, ", ", df_corr$CI_high)

val_dat_tri <- df_corr %>%
  rowwise() %>%
  mutate(pair = sort(c(Var1, Var2)) %>% paste(collapse = ",")) %>%
  group_by(pair) %>%
  distinct(pair, .keep_all = T)


val_dat_tri$Var1 <- factor(val_dat_tri$Var1, levels = unique(val_dat_tri$Var1))
val_dat_tri$Var2 <- factor(val_dat_tri$Var2, levels = unique(val_dat_tri$Var2))



p_sigma <- ggplot(val_dat_tri, aes(Var1, Var2)) +
      geom_tile(aes(fill = value), show.legend = TRUE) +
      geom_text(aes(label = CI), size = 5, nudge_y=-0.2) +
      scale_fill_gradient2(low = "#219ebc", high = "#ffb703", name = "Estimated Covariance") + 
      #scale_color_gradient2() + 
      geom_text(aes(label = round(value, 1)), size=6, nudge_y=0.1) +
      xlab("") + ylab("") +
      theme_minimal() + 
      theme(axis.text = element_text(size = 14, color = "black"))+
      scale_y_discrete(position = "right")

ggsave("./model2_corr.png", plot=p_sigma, width = 10, height = 7)


##POSTERIOR PREDICTIVE CHECKS##
sim_nets <- extract_network_draws(results$m_2, 1000)
dep_net <- results$m_2$data$dep_net


png("./plot.png")
simulated_network_checks(sim_nets, dep_net, "Triad_census")
dev.off()
simulated_network_checks(sim_nets, dep_net, "Indegree_distribution")
simulated_network_checks(sim_nets, dep_net, "Outdegree_distribution")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_baseline")
png("./plot.png")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_random")
dev.off()
