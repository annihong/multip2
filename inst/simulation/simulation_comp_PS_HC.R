
#library(multip2)
library(doSNOW)
# if (!check_simulation_requirement()){
#     message("make sure suggested packages are installed")
# }


#PARALLEL INFO
NUM_CORE = 16
TOTAL_ITER = 10
OUTPUT_PATH = "../simres/"
# RSTAN CONFIG:
CHAINS = 4
WARMUP = 500
THIN = 1
HYDRA = "M2"

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
cl <- parallel::makeCluster(NUM_CORE)  # Create a cluster with 4 workers
doSNOW::registerDoSNOW(cl)  # Register the cluster for use with foreach

#results <- foreach(iter=1:TOTAL_ITER, .packages="rstan") %dopar% {
results <- foreach(iter=1:TOTAL_ITER, .packages = "multip2") %dopar% {
    cat(paste0("we are on iteration ", iter, " hydra: ",HYDRA, "\n"), file = "simulation_log.txt", append=T)
    Sigma <- multip2::sample_Sigma(eta, ig_shape, ig_scale, t)
    sampled_params <- multip2::sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
    simulated_network <- multip2::simulate_network(n, t, sampled_params)


    M <- simulated_network #network to fit:

    if (sum(is.na(M))) {
        cat(paste0("The simulated network is NA...\n"), file = "simulation_log.txt", append=T)
    }


    m_fit <- Mp2Model(M)
    m_PS <- fit(m_fit, chains = CHAINS, warmup = WARMUP, thin = THIN, iter = WARMUP*2)
    m_HC <- fit(m_fit, chains = CHAINS, warmup = WARMUP, thin = THIN, iter = WARMUP*2, stan_file = "multiplex_p2_HC.stan")


    sim_result <- list(sampled_params=sampled_params, m_PS = m_PS, m_HC = m_HC)
    saveRDS(sim_result,file = paste0(OUTPUT_PATH,
                                 "_comp_res_", iter, "_out_of_", TOTAL_ITER,".Rds"))
    cat(paste0("Simulated result saved! ","\n"), file = "simulation_log.txt", append=T)
    rm(sim_result)
    rm(m_PS)
    rm(m_HC)
}

stopCluster(cl)
######## END FITTING THE MODEL IN RSTAN #########

######## ENDSIMULATIONS: #########
