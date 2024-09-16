file_path <- "/Users/annihong/Documents/Rprojects/fitted_models/"
files <- grep("*_study_2_models.Rds", list.files(file_path), value=TRUE)
class_ids <- sub("^([0-9]+).*", "\\1", files)

##LIBRARIES##
library(brms)
library(dplyr)
library(rstan)

##FUNCTIONS##
extract_model_info <- function(files, model_num, type = "fixed") {
    mean <- data.frame()
    se <- data.frame()
    rhat <- data.frame()
    for (file in files) {
        result <- readRDS(paste0(file_path, file))
        if (model_num == 1) {
            m_fit <- result$m_1
        } else {
            m_fit <- result$m_2
        }
        s <- multip2::summary.Mp2Model(m_fit)[[type]]
        mean <- rbind.data.frame(mean, s[,"mean"])
        se <- rbind.data.frame(se,s[,"se_mean"])
        rhat <- rbind.data.frame(rhat,s[,"Rhat"])
        colnames(mean) <- rownames(s)
        colnames(se) <- rownames(s)
        colnames(rhat) <- rownames(s)
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
model1Corr <- extract_model_info(files, 1, "correlation")
print(model1Corr$rhat)
model2Corr <- extract_model_info(files, 2, "correlation")
print(model2Corr$rhat)
model1sigma <- extract_model_info(files, 1, "sigma")
print(model1sigma$rhat)
model2sigma <- extract_model_info(files, 2, "sigma")
print(model2sigma$rhat)


## AGGREGATE RESULTS ##
m1fixed = aggregate_results(model1fixed, fixed_priors)
m2fixed = aggregate_results(model2fixed, fixed_priors)
m1Corr = aggregate_results(model1Corr, cov_priors)
m2Corr = aggregate_results(model2Corr, cov_priors)
m1sigma = aggregate_results(model1sigma, cov_priors)
m2sigma = aggregate_results(model2sigma, cov_priors)

model_results <- list(m1fixed = m1fixed, m2fixed = m2fixed, m1sigma = m1sigma, m2sigma=m2sigma, m1Corr = m1Corr, m2Corr = m2Corr)
save(model_results, file=paste0(file_path, "wave4_meta_results.RData"))

