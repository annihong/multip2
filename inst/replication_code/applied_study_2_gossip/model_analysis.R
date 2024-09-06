file_path <- "/home/annihong/projects/simres/fitted_models/"
files <- grep("*_study_2_models.Rds", list.files(file_path), value=TRUE)

##LIBRARIES##
library(brms)
library(dplyr)
library(rstan)

##FUNCTIONS##
extract_model_info <- function(files, model_num, pattern = "^PS_mu|^rho|^cross|fixed") {
    mean <- data.frame()
    se <- data.frame()
    rhat <- data.frame()
    for (file in files) {
        result <- readRDS(paste0(file_path, file))
        if (model_num == 1) {
            p2_fit <- result$m_1$fit_res$stan_fit
        } else {
            p2_fit <- result$m_2$fit_res$stan_fit
        }
        s <- rstan::summary(p2_fit)$summary
        res <- s[grep(pattern,rownames(s)),1:10]
        mean <- rbind.data.frame(mean, res[,1])
        se <- rbind.data.frame(se,res[,2])
        rhat <- rbind.data.frame(rhat,round(res[,10], 2))
        colnames(mean) <- rownames(res)
        colnames(se) <- rownames(res)
        colnames(rhat) <- rownames(res)
    }
     return(list(mean=mean, se=se, rhat=rhat))
}

aggregate_results <- function(model_res, priors, class = class_ids) {
    params <- colnames(model_res$mean)
    res <-list()
    for (i in 1:length(params)) {
        dat <- cbind.data.frame(mean = model_res$mean[,i], se = model_res$se[,i], class = class) 
        print(mean(dat$mean))
        m.brm <- brm(mean|se(se) ~ 1 + (1|class), data = dat, prior = priors, iter = 5000)
        res[[params[i]]] <- summary(m.brm)$fixed
    }
    res <- as.data.frame(do.call(rbind, res))
    return(res)
}

##META ANALYSIS##
fixed_priors <- c(prior(normal(0,20), class = Intercept),
            prior(cauchy(0,0.5), class = sd))

cov_priors <- c(prior(normal(0,100), class = Intercept),
            prior(cauchy(0,0.5), class = sd))

model1fixed <- extract_model_info(files, 1)
print(model1fixed$rhat)
model2fixed <- extract_model_info(files, 2)
print(model2fixed$rhat)
model1Sigma <- extract_model_info(files, 1, "^Corr")
print(model1Sigma$rhat)
model2Sigma <- extract_model_info(files, 2, "^Corr")
print(model2Sigma$rhat)


##LOAD IN THE STAN FIT OBJECTS##
fit1 <- results$m_1$fit_res$stan_fit
fit2 <- results$m_2$fit_res$stan_fit

##CONVERGENCE CHECK## NEEDS REVISION (figure out PS or HC, change the stan model used to estimation)

S <- ggs(fit2, family = "mu|rho|Corr|fixed")
S <- ggs(fit2, family = "^Corr\[(\d+),(?=\d+\]$)(\d+)\]$")
ggmcmc(S, file = paste0("diagnostics",".pdf"),  family = "mu|rho", param_page=6) #plot=c("traceplot", "running", "geweke"),
ptrace <- ggs_traceplot(S, family = "^Corr") +
      facet_wrap(~ Parameter, ncol = 6, scales = "free")
p <- ggs_grb(S, family = "^mu") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
ggs_geweke(S, family = "mu|rho") + facet_wrap(~ Parameter, ncol = 3, scales = "free")
ggs_running(S, family = "Sigma") + facet_wrap(~ Parameter, ncol = 6, scales = "free")
ggs_density(ggs(radon$s.radon, par_labels=P, family="sigma"))
ptrace <- ggs_traceplot(S, family = "mu")
ggs_running(S, family = "^mu")
ggs_Rhat(S, family = "^mu")
ggs_geweke(S, family = "^mu")
ggsave(ptrace)

##MODEL ESTIMATES##

summary.Mp2Model(results$m_1)
summary.Mp2Model(results$m_2)


##### RANDOM EFFECTS ###### NEED TO FIX THIS!

multi_result_random <- extract_model_estimates(fit2, params = c("Sigma"), CI = TRUE)
sig <- multi_result_random[, 2] * multi_result_random[, 3] > 0 
sig <- matrix(sig, ncol = 2 * t)
val <- matrix(multi_result_random[, 1], ncol = 2 * t)
upper <- matrix(multi_result_random[, 3], ncol = 2 * t)
lower <- matrix(multi_result_random[, 2], ncol = 2 * t)
vars <- c("pol_sender", "pol_receiver", "sci_sender", "sci_receiver", "influ_sender", "influ_receiver")

# Reshape the data
sig_dat <- reshape_matrix(sig, t, vars)
val_dat <- reshape_matrix(val, t, vars)
lower_dat <- reshape_matrix(lower, t, vars)
upper_dat <- reshape_matrix(upper, t, vars)
# Format credible interval values
val_dat$CI_low <- paste0("(", round(lower_dat$value, 1))
val_dat$CI_high <- paste0(round(upper_dat$value, 1), ")")
val_dat$CI <- paste0(val_dat$CI_low, ", ", val_dat$CI_high)

val_dat_tri <- val_dat %>%
  rowwise() %>%
  mutate(pair = sort(c(Var1, Var2)) %>% paste(collapse = ",")) %>%
  group_by(pair) %>%
  distinct(pair, .keep_all = T)

p_sigma <- ggplot(val_dat_tri, aes(Var1, Var2)) +
      geom_tile(aes(fill = value), show.legend = TRUE) +
      geom_text(aes(label = CI), size = 5, nudge_y=-0.2) +
      scale_fill_gradient2(low = "#219ebc", high = "#ffb703", name = "Estimated Covariance") + 
      #scale_color_gradient2() + 
      geom_text(aes(label = round(value, 1)), size=6, nudge_y=0.1) +
      xlab("") + ylab("") +
      theme_minimal() + 
      theme(axis.text = element_text(size = 14, color = "black"))

ggsave("model2_Sigma.png", plot=p_sigma, width = 10, height = 7)

##POSTERIOR PREDICTIVE CHECKS##
sim_nets <- extract_network_draws(results$m_2, 1000)
dep_net <- results$m_2$dep_net
simulated_network_checks(sim_nets, dep_net, "Triad_census")
simulated_network_checks(sim_nets, dep_net, "Indegree_distribution")
simulated_network_checks(sim_nets, dep_net, "Outdegree_distribution")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_baseline")
simulated_network_checks(sim_nets, dep_net, "multiplex_gof_random")
