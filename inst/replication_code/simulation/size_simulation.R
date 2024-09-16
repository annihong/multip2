
#library(multip2)
library(doSNOW)
# if (!check_simulation_requirement()){
#     message("make sure suggested packages are installed")
# }



OUTPUT_PATH = "/home/annihong/projects/simres/size_res/"
# RSTAN CONFIG:
CHAINS = 3
WARMUP = 1000
THIN = 1
# NETWORK INFO: 
n_seq = c(100, 150, 200)
t = 2
H = t*(t - 1)/2
#PARALLEL INFO
NUM_CORE = CHAINS*length(n_seq)



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

########  SIMULATING NETWORK FROM GROUND TRUTH #########

   
########  END SIMULATING NETWORK FROM GROUND TRUTH #########
cl <- parallel::makeCluster(NUM_CORE, outfile="LogSizeSim.txt")  # Create a cluster with 4 workers
doSNOW::registerDoSNOW(cl)  # Register the cluster for use with foreach

#results <- foreach(iter=1:TOTAL_ITER, .packages="rstan") %dopar% {
results <- foreach(iter=1:length(n_seq), .packages = "multip2") %dopar% {
    cat(paste0("we are on iteration ", iter, "\n"))
    n_seq = c(100, 150, 200)
    n = n_seq[iter]
    Sigma <- multip2::sample_Sigma(eta, ig_shape, ig_scale, t)
    sampled_params <- multip2::sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
    simulated_network <- multip2::simulate_network(n, t, sampled_params)


    M <- simulated_network #network to fit:

    if (sum(is.na(M))) {
        cat(paste0("The simulated network is NA...\n"))
    }


    m_fit <- multip2::Mp2Model(M)
    m_fit <-  multip2::fit(m_fit, chains = CHAINS, warmup = WARMUP, iter = WARMUP*2, network_sim = FALSE)
    cat(paste0("n = ", n))


    sim_result <- list(sampled_params=sampled_params,Mp2_fit = m_fit)
    saveRDS(sim_result,file = paste0(OUTPUT_PATH, "network_size_", n, ".Rds"))
    cat(paste0("Simulated result saved! ","\n"))
    rm(sim_result)
    rm(m_fit)
}

stopCluster(cl)
######## END FITTING THE MODEL IN RSTAN #########

######## ENDSIMULATIONS: #########
