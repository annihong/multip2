
#library(multip2)
library(doSNOW)
# if (!check_simulation_requirement()){
#     message("make sure suggested packages are installed")
# }
path <- "/home/annihong/projects/simres/"
simulations <- list.files(paste0(path, "sim_res/new/"))
done_ids <- na.omit(stringr::str_extract(simulations, "^[0-9]+"))
iter_ids <- setdiff(1:1000, done_ids)
#simulations <- grep("1000.Rds$", simulations, value = TRUE)[1:2]
#simulations <- grep("test_", list.files(paste0(path, "sim_res/")), value=TRUE)[1:20]
prior_res <- readRDS("/home/annihong/projects/simres/analysis_res/prior_res1000.rds")
OUTPUT_PATH = "/home/annihong/projects/simres/sim_res/new/"
# iter_ids = 901:1000
#PARALLEL INFO
NUM_CORE = 20

# RSTAN CONFIG:
CHAINS = 2
WARMUP = 1000
THIN = 1
EXPERIMENT = "sim:rest"


# NETWORK INFO: 
n = 30
t = 2
H = t*(t - 1)/2


########  END SIMULATING NETWORK FROM GROUND TRUTH #########
cl <- parallel::makeCluster(NUM_CORE, outfile=paste0(OUTPUT_PATH, "log/", EXPERIMENT, "LogSim.txt"))  # Create a cluster with 4 workers
doSNOW::registerDoSNOW(cl)  # Register the cluster for use with foreach

#results <- foreach(iter=1:TOTAL_ITER, .packages="rstan") %dopar% {
results <- foreach(i=1:length(iter_ids), .packages = "multip2") %dopar% {
    iter = iter_ids[i]
    sampled_params <- prior_res[[iter]]$sampled_params
    sampled_params$iter = iter
    M <- prior_res[[iter]]$simulated_network #network to fit:
    m_fit <- multip2::Mp2Model(M)
    m_fit <-  multip2::fit(m_fit, chains = CHAINS, warmup = WARMUP, thin = THIN, iter = WARMUP*2, network_sim = TRUE, seed = iter, par = "x_beta", include = FALSE, stan_file = "multiplex_p2_revert.stan")
    # m_fit <-  multip2::fit(m_fit, chains = 1, warmup =10, thin = 1, iter = 1000, network_sim = FALSE, prior_sim = TRUE, stan_file = "multiplex_p2_low_mem_prior_sim_test.stan")


    sim_result <- list(sampled_params=sampled_params,Mp2_fit = m_fit)
    saveRDS(sim_result,file = paste0(OUTPUT_PATH, iter, "_", EXPERIMENT, "draws = ", CHAINS * WARMUP, ".Rds"))
    rm(sim_result)
    rm(m_fit)
}

stopCluster(cl)
######## END FITTING THE MODEL IN RSTAN #########

######## ENDSIMULATIONS: #########
