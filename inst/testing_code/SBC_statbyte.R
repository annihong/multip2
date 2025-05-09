devtools::install_github("hyunjimoon/SBC")
library(SBC)


####### Generator #########
# A generator function should return a named list containing elements "variables" and "generated"
statbyte_example_generator_single <- function(n, K, mu_0){  # sim_N is the number of data points we are generating
   N = n*K
   group_idx = rep(1:K, each = n)
   #sampling from the prior
   mu = rnorm(1, mu_0, 5)
   sigma = extraDistr::rinvgamma(1, 1,1)
   beta_ks = rnorm(K, mu, sigma)
    mu_ps = mu + mean(beta_ks)  

    #simulate data from the prior draw
    prob_logit = mu + beta_ks[group_idx]
    prob = exp(prob_logit)/(1 + exp(prob_logit))
    y = rbinom(N, 1, prob)

    stan_data = list(
        N = N,
        K = K,
        y = y,
        group_idx = group_idx,
        sample_prior = 0
    )

  list(
    variables = list(
        mu = mu,
        sigma = sigma,
        mu_ps = mu_ps
    ),
    generated = stan_data
  )
}


set.seed(54882235)
n_sims <- 100  # Number of SBC iterations to run

statbyte_example_generator <- SBC_generator_function(statbyte_example_generator_single, n = 30, K = 5, mu_0 = 0)


dataset <- generate_datasets(
  statbyte_example_generator , 
  n_sims)


rstan_model <- rstan::stan_model("/home/annihong/projects/multip2/inst/stan/hierarchical_logit.stan")


backend <- SBC_backend_rstan_sample( rstan_model, iter = 2000, warmup = 1000, chains = 2)  

results <- compute_SBC(dataset, backend, 
                    cache_mode = "results", 
                    cache_location = file.path("/home/annihong/projects/simres/multilevel_p2", "results"))

# results <- readRDS("/home/annihong/projects/simres/multilevel_p2/results.rds")
results$stats
p1 <- plot_rank_hist(results, bins=20)
p2 <- plot_ecdf(results)

p3 <- plot_coverage(results)
p4 <- plot_coverage_diff(results)
p5 <- plot_ecdf_diff(results)


p <- gridExtra::grid.arrange(p1, p2, p3, p4,p5, nrow=5)
ggplot2::ggsave("/home/annihong/projects/multip2/inst/testing_code/SBC_statbyte.png", p, width = 10, height = 15, dpi = 300)







