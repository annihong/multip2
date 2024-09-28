
#library(multip2)
library(doSNOW)
# if (!check_simulation_requirement()){
#     message("make sure suggested packages are installed")
# }

prior_res <- readRDS("/home/annihong/projects/simres/analysis_res/prior_res.rds")
rho_s <- sapply(prior_res, function(x) x$sampled_params$rho)
mu_s <- sapply(prior_res, function(x) x$sampled_params$mu)
#PARALLEL INFO
NUM_CORE = 32
TOTAL_ITER = length(prior_res)
OUTPUT_PATH = "/home/annihong/projects/simres/sim_res/"
# RSTAN CONFIG:
CHAINS = 1
WARMUP = 800
THIN = 1
HYDRA = "test"
CURRENT = 1
END_ITER = CURRENT + TOTAL_ITER - 1
# NETWORK INFO: 
n = 30
t = 2
H = t*(t - 1)/2


# true mean of the priors
#Sigma = diag(2*t)
mu_0 = rep(0,t)
rho_0 = rep(0,t)
cross_mu_0 = rep(0,H)
cross_rho_0 = rep(0,H)
eta = 2
Sigma_Dist <- "LJK" 
ig_shape = 3
ig_scale = 50

TIME = Sys.time()

cat(as.character(TIME) , file = "simulation_log.txt", append=T)
cat("\n" , file = "simulation_log.txt", append=T)
cat(paste0("There are a total of ", TOTAL_ITER, " simulations, and # iterations = ",2*WARMUP," Let's hope for the best!", "hydra: ",HYDRA, "\n"), file = "simulation_log.txt", append=T)
cat(paste0("prior info: eta = ", eta, ", shape = ", ig_shape, ", scale= ",ig_scale, "\n"), file = "simulation_log.txt", append=T)

######## END USER INPUT #########

######## BEGIN SIMULATING NETWORK FROM GROUND TRUTH #########



########  SIMULATING NETWORK FROM GROUND TRUTH #########

   
########  END SIMULATING NETWORK FROM GROUND TRUTH #########
cl <- parallel::makeCluster(NUM_CORE, outfile="LogSim.txt")  # Create a cluster with 4 workers
doSNOW::registerDoSNOW(cl)  # Register the cluster for use with foreach

#results <- foreach(iter=1:TOTAL_ITER, .packages="rstan") %dopar% {
results <- foreach(iter=CURRENT:END_ITER, .packages = "multip2") %dopar% {
    cat(paste0("we are on iteration ", iter, " hydra: ",HYDRA, "\n"), file = "simulation_log.txt", append=T)
    # Sigma <- multip2::sample_Sigma(eta, ig_shape, ig_scale, t)
    # sampled_params <- multip2::sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
    # simulated_network <- multip2::simulate_network(n, t, sampled_params)

    sampled_params <- prior_res[[iter]]$sampled_params
    M <- prior_res[[iter]]$simulated_network #network to fit:

    if (sum(is.na(M))) {
        cat(paste0("The simulated network is NA...\n"), file = "simulation_log.txt", append=T)
    }


    m_fit <- multip2::Mp2Model(M)
    m_fit <-  multip2::fit(m_fit, chains = CHAINS, warmup = WARMUP, thin = THIN, iter = WARMUP*2, network_sim = FALSE)
    # m_fit <-  multip2::fit(m_fit, chains = 1, warmup =10, thin = 1, iter = 1000, network_sim = FALSE, prior_sim = TRUE, stan_file = "multiplex_p2_low_mem_prior_sim_test.stan")


    sim_result <- list(sampled_params=sampled_params,Mp2_fit = m_fit)
    saveRDS(sim_result,file = paste0(OUTPUT_PATH, HYDRA,
                                 "_", iter, "_out_of_", END_ITER, "draws = ", CHAINS * WARMUP, ".Rds"))
    cat(paste0("Simulated result saved! ","\n"), file = "simulation_log.txt", append=T)
    rm(sim_result)
    rm(m_fit)
}

stopCluster(cl)
######## END FITTING THE MODEL IN RSTAN #########

######## ENDSIMULATIONS: #########
