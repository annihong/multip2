######## SAMPLE FROM PRIOR GIVEN GROUND TRUTH MEANS #########
# #' Generate a covariance matrix using the LKJ decomposition method
# #'
# #' @param None
# #' 
# #' @return A covariance matrix
# Sigma_cor<- function(){
#     sigma <- sigma_folded_normal(0,10,2)
#     #sigma = sigma_cauchy(1,2)
#     L_corr <- LKJ_decomp(2, 2)
#     Sigma = (diag(sigma) %*% L_corr) %*% t(diag(sigma) %*% L_corr) 
#     return(Sigma)
# }

#' Generate a correlation matrix using the LKJ decomposition method
#'
#' @param eta Parameter for the LKJ distribution
#' @param t number if networks 
#' 
#' @return A correlation matrix
LKJ_decomp <- function(eta, t){

    my_lkjcorr_fun = "
    functions {
    // generate lkj correlation matrix (R)
    matrix my_lkj_corr_rng(int K, real eta) {
        return lkj_corr_rng(K, eta);
    }
    
    // generate cholesky factor L_corr of a correlation matrix R
    matrix my_lkj_corr_chol_rng(int K, real eta){
    return lkj_corr_cholesky_rng(K, eta);
    }
    
    // perform triangular matrix multiplication L*L^T
    matrix my_multiply_lower_tri_self_transpose(matrix L){
        return multiply_lower_tri_self_transpose(L);
    }
    }
    "
    rstan::expose_stan_functions(rstan::stanc(model_code = my_lkjcorr_fun))
    L_corr = my_lkj_corr_chol_rng(K = 2*t, eta = eta)
    return(L_corr)
}

#' Generate a vector of random numbers from a folded normal distribution
#'
#' @param mean Mean of the folded normal distribution
#' @param sd Standard deviation of the folded normal distribution
#' @param t Number of networks
#' 
#' @return A vector of random numbers
sigma_folded_normal <- function(mean, sd, t){
    return(greybox::rfnorm(2*t, mean, sd))
}

#' Generate a vector of random numbers from an inverse gamma distribution
#'
#' @param shape Shape parameter of the inverse gamma distribution
#' @param scale Scale parameter of the inverse gamma distribution
#' @param t Number of networks
#' 
#' @return A vector of random numbers
sigma_inver_gamma <- function(shape, scale, t) {
    return(extraDistr::rinvgamma(2*t, shape, scale))
}

#' Generate a vector of random numbers from a Cauchy distribution
#'
#' @param k Scale parameter of the Cauchy distribution
#' @param t Number of networks
#' 
#' @return A vector of random numbers
sigma_cauchy <- function(k, t){
    sigma_unif = runif(2*t, 0, pi/2)
    sigma = k * tan(sigma_unif)
    return(sigma)
}

#' Generate a covariance matrix using the LKJ decomposition method
#'
#' @param sigma A vector of standard deviations
#' @param L_corr A correlation matrix
#' 
#' @return A covariance matrix
Sigma_LJK <- function(sigma, L_corr) {
    Sigma = (diag(sigma) %*% L_corr) %*% t(diag(sigma) %*% L_corr)
    return(Sigma) 
    
}
# Sigma_inwish <- function(t, nu=2*t + 1, S=diag(2*t) ) {
#     Sigma <- LaplacesDemon::rinvwishart(nu, S)
#     return(Sigma)
# }


#given a list of functions that corresponds to a prior distribution
#' sample_prior function
#'
#' This function generates samples from the priors.
#'
#' @param n The number of samples to generate.
#' @param t The number of networks.
#' @param mu_0 The prior mean for mu.
#' @param rho_0 The prior mean for rho.
#' @param cross_mu_0 The prior mean for cross_mu.
#' @param cross_rho_0 The prior mean for cross_rho.
#' @param Sigma The covariance matrix.
#'
#' @return A list containing the generated samples for mu, rho, cross_mu, cross_rho, C, and Sigma.

sample_prior <- function(n, t, mu_0, rho_0, cross_mu_0, cross_rho_0, Sigma){
    #simulate from the priors 
    mu = sapply(mu_0, rnorm, n=1, sd=10)
    rho = sapply(rho_0, rnorm, n=1, sd=10)
    cross_mu = sapply(cross_mu_0, rnorm, n=1, sd=10)
    cross_rho = sapply(cross_rho_0, rnorm, n=1, sd=10)
    
    C = mvtnorm::rmvnorm(n=n,mean=rep(0, 2*t),sigma=Sigma)
    
    return(list(mu=mu, rho=rho, cross_mu=cross_mu, cross_rho=cross_rho, C=C,
                 Sigma=Sigma))
}

#' @export
sample_Sigma <- function(eta, ig_shape, ig_scale, t) {
#   if (Sigma_Dist == "LJK") {
#       L_corr = LKJ_decomp(eta,t)
#       #sigma = sigma_folded_normal(sigma_0,sigma_sd,t)
#       sigma = sigma_inver_gamma(ig_shape, ig_scale, t)
#       Sigma = Sigma_LJK(sigma, L_corr)
#   } else {
#       Sigma = Sigma_inwish(t, nu)
#   }
    L_corr = LKJ_decomp(eta,t)
    sigma = sigma_inver_gamma(ig_shape, ig_scale, t)
    Sigma = Sigma_LJK(sigma, L_corr)
  return(Sigma)
}


######## END SAMPLE FROM PRIOR GIVEN GROUND TRUTH MEANS #########

######## HELPER FUNCTIONS FOR SIMULATE NETWORKS GIVEN PARAMETERS #########

#' Given number of actors, total number of networks, and the var-cov matrix
#' return a list of sender (alpha) \in R^{n x T} and receiver effects (beta)
#'
#' @param n_actor Number of actors
#' @param T_network Total number of networks
#' @param Sigma Variance-covariance matrix
#' @param supplied_C Boolean indicating whether the var-cov matrix is supplied or generated
#' @return A list containing the sender effects (alpha) and receiver effects (beta)
#' @export
sim_random_effects <- function(n_actor, T_network, Sigma, supplied_C=FALSE) {
    if (supplied_C) {
        C = Sigma
    } else {
        C = mvtnorm::rmvnorm(n_actor, mean = rep(0, nrow(Sigma)), sigma = Sigma)
    }
    alpha <- as.matrix(C[,seq(1, 2 * T_network, 2)])
    beta <- as.matrix(C[,seq(2, 2 * T_network, 2)])
    return(list(alpha = alpha, beta = beta))
}


calc_within_term <- function(i,j,k,T_network,alpha,beta,mu,rho,m1_all,m2_all){
    calc_within_term_t <- function(i,j,k,t){
        m1 = m1_all[k,t]
        m2 = m2_all[k,t]
        within_terms = m1*(alpha[i,t] + beta[j,t] + mu[t])
        within_terms = within_terms + m2*(alpha[j,t] + beta[i,t]+ mu[t])
        within_terms = within_terms + m1*m2*(rho[t])
        return(within_terms)
    }
    within_terms = 0
    for (l in 1:T_network) { #adding up contribution from all networks T
        within_terms = within_terms + calc_within_term_t(i,j,k,l)
    }
    return(within_terms)
}


calc_cross_term <- function(i,j,k,H,cross_mu,cross_rho, network_pairs, m1_all, m2_all){
    calc_cross_term_h <- function(i,j,k,h){
        t_a = network_pairs[[h]][1]
        t_b = network_pairs[[h]][2]

        m1_a = m1_all[k,t_a]
        m2_a = m2_all[k,t_a]
        m1_b = m1_all[k,t_b]
        m2_b = m2_all[k,t_b]

        cross_terms = (m1_a * m1_b + m2_a * m2_b) * cross_mu[h]
        cross_terms = cross_terms + (m1_a * m2_b + m2_a * m1_b) * cross_rho[h]
    }

    cross_terms = 0

    for (h in 1:H) { #adding up contribution from all H pairs of networks
        cross_terms = cross_terms + calc_cross_term_h(i,j,k,h)
    }
    return(cross_terms)
}

calc_m1 <- function(score) {
    m1 = (score == 2 | score == 4)
    return(m1)
}

calc_m2 <- function(score) {
    m2 = (score == 3 | score == 4)
    return(m2)
}

check_m1m2_calculation <- function(){
    for (k in outcomes) {
        for (j in 1:t) {
            score = one_d_res[k,j]
            m1 = (score == 2 | score == 4)
            m2 = (score == 3 | score == 4)
            assertthat::assert_that(are_equal(m1, m1_all[k,j]))
            assertthat::assert_that(are_equal(m2, m2_all[k,j]))
        }
    }
}

check_true_param_dim <- function() {
    #assert_that(are_equal(dim(Sigma), c(2*t, 2*t)))
    assertthat::assert_that(are_equal(length(mu), t))
    assertthat::assert_that(are_equal(length(rho), t))
    assertthat::assert_that(are_equal(length(cross_mu), H))
    assertthat::assert_that(are_equal(length(cross_rho), H))
}



######## END HELPER FUNCTIONS FOR SIMULATE NETWORKS GIVEN PARAMETERS #########


######## SIMULATE NETWORKS GIVEN PARAMETERS #########
calculate_probs_dyad_outcomes <- function(n, t, params){
    # network stuff setup
    H = t*(t - 1)/2
    N = n*(n - 1)/2 # number of dyads
    K = 2^(2*t)
    outcomes = 1:K
    one_d_res = y_nd_ij_to_y_1d_ij(outcomes, t)

    dyads = get_dyads(n)
    network_pairs = get_dyads(t)

    m1_all = apply(one_d_res, 1:2, calc_m1)
    m2_all = apply(one_d_res, 1:2, calc_m2) 

    # calculate probability for each outcome on all the dyads
    sim_prob_list = data.frame()
    for (dyad in dyads) {
        i <- dyad[1]
        j <- dyad[2]

        res = numeric(length(outcomes))

        for (k in outcomes) {
            within_terms <-
            calc_within_term(i,j,k,t, params$alpha, params$beta, params$mu, params$rho, m1_all, m2_all)
            if (t < 2) {
                cross_terms <- 0
            } else {
                cross_terms <- calc_cross_term(i,j,k,H,params$cross_mu, params$cross_rho, network_pairs, m1_all, m2_all)
                res[[k]] <- within_terms + cross_terms
            }
            
        }
        res <- exp(res)
        sim_prob <- res/(sum(res))
        sim_prob_list <- rbind(sim_prob_list, sim_prob)
    }

    colnames(sim_prob_list) <- outcomes
    return(sim_prob_list)
}

#simulate one multiplex network (with t layers) given parameters
#' @export
simulate_network <- function(n, t, params){
    #helper for sampleing from the prob list 
    K <- 2^(2*t)
    outcomes <- 1:K

    sample_y <- function(probs) {
        #print(probs)
        #print(probs)
        sample(outcomes, 1, replace = FALSE, prob = probs)
    }

    #obtain actor effects
    actor_eff <- sim_random_effects(n, t, params$C, TRUE) #use C instead of Sigma
    params$alpha <- actor_eff$alpha
    params$beta <- actor_eff$beta

    probs <-  calculate_probs_dyad_outcomes(n, t, params)
    
    #print(probs)
    sim_ys <- list()
    for (i in 1:nrow(probs)){
        sim_ys[[i]] <- tryCatch(sample_y(probs[i,]),
        error = function(e)
        paste("bad prior values")
            ) 
    }
    if ("bad prior values" %in% sim_ys){
        networks = NA
    } else {
        networks <- dyads_to_matrix_nd(sim_ys,t,n)
    }
    #print(networks)
    return(networks)
}

######## END SIMULATE NETWORKS GIVEN PARAMETERS #########



# descriptive stats based goodness of fit
# post processing the generated results
# M_net <- lapply(M, igraph::graph_from_adjacency_matrix)
# descriptive_statistics <- descriptive_stats(M_net)

# observed_stats <- descriptive_statistics
# observed_stats_df <- data.frame(var=names(observed_stats), sim_stats= observed_stats)


# sim_sample_size = 100
# sim_nws <- tail(rstan::extract(p2_fit)$y_tilde, sim_sample_size)
# sim_nws_igraphs <- dyads_to_matrix_list(sim_nws, n, t, "igraph") #convert the simulated dyad data to igraph objects 

# #calculate the descriptive stats for all the simulated data 
# stats <- descriptive_stats_list(sim_nws_igraphs, avg=F)
# basic_stats <- gather(data.frame(stats[,1:6]), key="var", value="sim_stats")
# varcov_stats <- gather(data.frame(stats[,7:22]), key="var", value="sim_stats")


# sigma_seq <- replicate(1000, sigma_cauchy(k=10, t=2))
# apply(sigma_seq, 1, mean)
# apply(sigma_seq, 1, sd)
# sigma_seq <- replicate(1000, sigma_normal(10,5,2))
# apply(sigma_seq, 1, mean)
# apply(sigma_seq, 1, sd)

# apply(replicate(10, Sigma_cor()), 1:2, mean)
# apply(replicate(100, Sigma_cor()), 1:2, sd)

# L = nimble::rlkj_corr_cholesky(n = 1, eta=300, p=2*t)
# 
# #' sim_L function
# #'
# #' This function generates a correlation matrix using the LKJ correlation distribution.
# #'
# #' @param eta The shape parameter of the LKJ correlation distribution.
# #' @param t The number of networks.
# #' @importFrom nimble rlkj_corr_cholesky
# #' @return A correlation matrix generated using the LKJ correlation distribution.
# sim_L <- function(eta,t) {
#     L = nimble::rlkj_corr_cholesky(n = 1, eta=eta, p=2*t)
#     cor <- stats::cov2cor(L %*% t(L))
#     return(cor)
# }

# sim_L_corr <- function() {
#     # Set parameters
#     t <- 2  # Degree of correlation among each group of 2 variables
#     n = 100
#     eta_seq <- c(0.5, 1, 2, 5, 10, 30, 50)  # Sequence of eta values to simulate

#     # Generate list of averaged correlation matrices for each eta value
#     avg_cor_matrices <- list()
#     sd_cor_matrices <- list()
#     for (i in 1:length(eta_seq)) {
#     cor_matrices <- replicate(n, sim_L(eta_seq[i], t), simplify = FALSE)
#     avg_cor_matrix <- Reduce("+", cor_matrices) / n
#     avg_cor_matrices[[i]] <- round(avg_cor_matrix,2)
#     sd_cor_matrix <- apply(simplify2array(cor_matrices), MARGIN = c(1, 2), sd)
#     sd_cor_matrices[[i]] <- sd_cor_matrix  
#     }

#     names(avg_cor_matrices) <- eta_seq
#     names(sd_cor_matrices) <- eta_seq
#     avg_cor_matrices
#     sd_cor_matrices
# }

#' Check simulation requirements
#'
#' This function checks if the required packages for running simulations are loaded.
#' It checks for the presence of the "doSNOW", "greybox", "mvtnorm", and "extraDistr" packages.
#'
#' @import doSNOW
#' @return A logical value indicating whether all the required packages are loaded.
check_simulation_requirement <- function() {
    return(requireNamespace("doSNOW", quietly = TRUE)
            & requireNamespace("greybox", quietly = TRUE)
            & requireNamespace("mvtnorm", quietly = TRUE) 
            & requireNamespace("extraDistr", quietly = TRUE)
            )
}