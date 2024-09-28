
#library(multip2)
library(doSNOW)
# if (!check_simulation_requirement()){
#     message("make sure suggested packages are installed")
# }


#PARALLEL INFO
NUM_CORE = 32
TOTAL_ITER = 10
OUTPUT_PATH = "/home/annihong/projects/simres/sim_res/"
# RSTAN CONFIG:
CHAINS = 2
WARMUP = 500
THIN = 1
HYDRA = "G2"
CURRENT = 1000
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

check_na <- function(M) {
    has_na <- is.na(sum(sapply(M, sum)))
    return(has_na)
}

get_sample <- function() {
    Sigma <- multip2::sample_Sigma(eta, ig_shape, ig_scale, t)
    res <- multip2::sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
    res$simulated_network <- multip2::simulate_network(n, t, sampled_params)
    res$is_na <- check_na(res$simulated_network)
    return(res)
}



res <- replicate(200, get_sample())
mu_s <- sapply(1:200, function(i) res[,i]$mu)
rho_s <- sapply(1:200, function(i) res[,i]$rho)



Sigma <- multip2::sample_Sigma(eta, ig_shape, ig_scale, t)
sampled_params <- multip2::sample_prior(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma)
return(sampled_params$rho)
simulated_network <- multip2::simulate_network(n, t, sampled_params)