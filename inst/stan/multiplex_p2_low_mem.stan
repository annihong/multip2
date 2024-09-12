functions  {
    int[] find_start_end(int[] param_count, int t){
    
        int res[2];
        res[1] = 1;
        res[2] = 0;
        
      if (size(param_count) < t) {
        reject("selected layer index (for pair index) larger than the total number of layers; found layer_t (h) =", t);
      }
        //for starting index
        for (i in 1:(t-1)) {
          res[1] = res[1] + param_count[i];
        }
        
        //for ending index
        for (i in 1:t) {
          res[2] = res[2] + param_count[i];
        }
        return res;
      }
  }

data {
  int<lower=0> prior_sim; //boolean for whether to simulate from priors (1) or fit the actual model (0)
  int<lower=0> n; //number of actors
  int<lower=0> N; // number of total dyads
  int<lower=0> T; // number of networks
  int<lower=0> H; // number of pairs of layers of networks
  int<lower=0> layer_pairs[H != 0 ? H : 1, 2]; //a list that stores all the network indices of the pairs of layers
  // if H > 0 (at least two networks), then H rows in the layer pairs else 1 row
  int<lower=0> N_obs;
  int<lower=1, upper=N> obs_idx[N_obs];


    //outcomes
  int<lower=0> K; //number of outcomes on a dyad := 2^2T
  int<lower=1, upper=K> y_obs[N_obs]; //y[N] \in {1,...,2^{2T}}^N, actual observed outcomes

  
  //adding covariates to multiplex p2
  int<lower=0> D_within[T, 4]; // X[t, 1:4]  of covariates for mu_t, rho_t, alpha_t, beta_t
  int<lower=0> D_cross[H, 2]; // X[h, 1:2] of covar for cross_mu_h, cross_rho_h
  
  vector[sum(D_within[,1])] mu_covariates[n, n]; // for mu  
  vector[sum(D_within[,2])] rho_covariates[n, n]; // for rho
  matrix[n,sum(D_within[,3])] alpha_covariates; // for alpha
  matrix[n,sum(D_within[,4])] beta_covariates; // for beta
  
  vector[H != 0 ? sum(D_cross[,1]) : 0] cross_mu_covariates[n, n];
  vector[H != 0 ? sum(D_cross[,2]) : 0] cross_rho_covariates[n, n];


  // priors
  vector[T] mu_mean_prior;
  vector[T] mu_sd_prior;
  vector[T] rho_mean_prior;
  vector[T] rho_sd_prior;

  vector[H] cross_mu_mean_prior;
  vector[H] cross_mu_sd_prior;
  vector[H] cross_rho_mean_prior;
  vector[H] cross_rho_sd_prior;

  real scale_alpha_prior;
  real scale_beta_prior;
  real LJK_eta_prior;

  // priors for the covariates

  vector[sum(D_within[,1])] mu_covariates_sd_prior;
  vector[sum(D_within[,1])] mu_covariates_mean_prior;

  vector[sum(D_within[,2])] rho_covariates_sd_prior;
  vector[sum(D_within[,2])] rho_covariates_mean_prior;

  vector[sum(D_within[,3])] alpha_covariates_sd_prior;
  vector[sum(D_within[,3])] alpha_covariates_mean_prior;

  vector[sum(D_within[,4])] beta_covariates_sd_prior;
  vector[sum(D_within[,4])] beta_covariates_mean_prior;

  vector[H != 0 ? sum(D_cross[,1]) : 0] cross_mu_covariates_sd_prior;
  vector[H != 0 ? sum(D_cross[,1]) : 0] cross_mu_covariates_mean_prior;

  vector[H != 0 ? sum(D_cross[,2]) : 0] cross_rho_covariates_sd_prior;
  vector[H != 0 ? sum(D_cross[,2]) : 0] cross_rho_covariates_mean_prior;
}

parameters {
  vector[T] mu; // vector of size T for within-network density
  vector[T] rho; // vector of size T for within-network reciprocity
  vector[H] cross_mu; // size H storing cross-network density for each pair of network
  vector[H] cross_rho; // size H storing cross-network reciprocity for each pair of network
  
  //fixed effects:
  vector[sum(D_within[,1])] mu_fixed_coef;  
  vector[sum(D_within[,2])] rho_fixed_coef; 
  vector[sum(D_within[,3])] alpha_fixed_coef; 
  vector[sum(D_within[,4])] beta_fixed_coef; 
  vector[H != 0 ? sum(D_cross[,1]) : 0] cross_mu_fixed_coef; 
  vector[H != 0 ? sum(D_cross[,2]) : 0] cross_rho_fixed_coef;
  
  //random effects:
  //cov_matrix[2*T] Sigma; // cov matrix used to draw the random actor effects
  vector<lower=0>[2*T] sigma;
  cholesky_factor_corr[2*T] L_corr;
  matrix[2*T, n] z;
  
}

transformed parameters{
  
}

model {
  matrix[N,K] x_beta;
  
  // start of the x_beta calculation
  
  {
    
    matrix[n,2*T] C; // for each actor, there are 2 * T number of random actor effects (two per network)
    matrix[n,T] alpha;
    matrix[n,T] beta;
    real mu_ij = 0;
    real mu_ji = 0;
    real rho_ij = 0;
    real cross_mu_ij = 0;
    real cross_rho_ij = 0;
    vector[n] alpha_fixed = rep_vector(0,n);
    vector[n] beta_fixed = rep_vector(0,n);
  
    C = (diag_pre_multiply(sigma, L_corr) * z)';
    for (t in 1:T){
      int idx_a[2] = find_start_end(D_within[,3],t);
      int idx_b[2] = find_start_end(D_within[,4],t);


      if (D_within[t,3] > 0) {
        alpha_fixed = alpha_covariates[,idx_a[1]:idx_a[2]] * alpha_fixed_coef[idx_a[1]:idx_a[2]];
      }
      
      if (D_within[t,4] > 0) {
        beta_fixed = beta_covariates[,idx_b[1]:idx_b[2]] * beta_fixed_coef[idx_b[1]:idx_b[2]];
      }
      
      alpha[,t] = C[,1 + 2 * (t - 1)] + alpha_fixed;
      beta[,t] = C[,2 + 2 *(t - 1)] + beta_fixed;
    }
    
    {int counter = 1;
    matrix[T,2] M; 
      for (i in 1:n) {
        for (j in i:n) {
          if (i == j) {
            continue;
          }
          if (counter )
          for (k in 1:K){
  
            real within_terms = 0;
            real cross_terms = 0;
            
            for (t in 1:T) {
              real nt = ceil(k / 4^(t - 1));
              real score = fmod(nt,4);
              int idx_mu[2] = find_start_end(D_within[,1],t);
              int idx_rho[2] = find_start_end(D_within[,2],t);

              M[t,] = [(score == 2 || score == 0),(score == 3 || score == 0)];
              
              mu_ij = D_within[t,1] > 0 ? dot_product(mu_covariates[i,j][idx_mu[1]:idx_mu[2]], mu_fixed_coef[idx_mu[1]:idx_mu[2]]) : 0;
              mu_ji = D_within[t,1] > 0 ? dot_product(mu_covariates[j,i][idx_mu[1]:idx_mu[2]], mu_fixed_coef[idx_mu[1]:idx_mu[2]]) : 0;
              rho_ij = D_within[t,2] > 0 ? dot_product(rho_covariates[i,j][idx_rho[1]:idx_rho[2]], rho_fixed_coef[idx_rho[1]:idx_rho[2]]) : 0;
            

            //print("rho dot:", dot_product(rho_covariates[i,j][idx_rho[1]:idx_rho[2]], rho_fixed_coef[idx_rho[1]:idx_rho[2]]));
            
              within_terms += M[t,1]*(alpha[i,t] + beta[j,t] + mu[t] + mu_ij); 
              within_terms += M[t,2]*(alpha[j,t] + beta[i,t] + mu[t] + mu_ji); 
              within_terms += M[t,1]*M[t,2]*(rho[t] + rho_ij);
              //print("reciprocated: ", M[t,1]*M[t,2]); 

              
            }
            
            for (h in 1:H){
              int net_a = layer_pairs[h,1];
              int net_b = layer_pairs[h,2];
              
              int idx_mu[2] = find_start_end(D_cross[,1],h);
              int idx_rho[2] = find_start_end(D_cross[,2],h);


              cross_mu_ij = D_cross[h,1] > 0 ? dot_product(cross_mu_covariates[i,j][idx_mu[1]:idx_mu[2]],cross_mu_fixed_coef[idx_mu[1]:idx_mu[2]]): 0;
              cross_rho_ij = D_cross[h,2] > 0 ? dot_product(cross_rho_covariates[i,j][idx_rho[1]:idx_rho[2]], cross_rho_fixed_coef[idx_rho[1]:idx_rho[2]]) : 0;
                
                
              cross_terms += (M[net_a,1] * M[net_b,1] + M[net_a,2] * M[net_b,2]) * (cross_mu[h] + cross_mu_ij);
              cross_terms += (M[net_a,1] * M[net_b,2] + M[net_a,2] * M[net_b,1]) * (cross_rho[h] + cross_rho_ij);
            }
            
            x_beta[counter,k] = within_terms + cross_terms;
           
          }
          counter += 1;
        }
      }
    }
  }
  // end of the x_beta calculation


  for (i in 1:T) {
    mu[i] ~ normal(mu_mean_prior[i],mu_sd_prior[i]);
    rho[i] ~ normal(rho_mean_prior[i],rho_sd_prior[i]);
  }

  for (i in 1:H) {
    cross_mu[i] ~ normal(cross_mu_mean_prior[i],cross_mu_sd_prior[i]);
    cross_rho[i] ~ normal(cross_rho_mean_prior[i],cross_rho_sd_prior[i]);
  }
  to_vector(z) ~ std_normal();
  L_corr ~ lkj_corr_cholesky(LJK_eta_prior);
  sigma ~ inv_gamma(scale_alpha_prior,scale_beta_prior);


  for (i in 1:sum(D_within[,1])) {
    mu_fixed_coef[i] ~ normal(mu_covariates_mean_prior[i],mu_covariates_sd_prior[i]);
  }
  
  for (i in 1:sum(D_within[,2])) {
    rho_fixed_coef[i] ~ normal(rho_covariates_mean_prior[i],rho_covariates_sd_prior[i]);
  }
  
  for (i in 1:sum(D_within[,3])) {
    alpha_fixed_coef[i] ~ normal(alpha_covariates_mean_prior[i],alpha_covariates_sd_prior[i]);
  }
  
  for (i in 1:sum(D_within[,4])) {
    beta_fixed_coef[i] ~ normal(beta_covariates_mean_prior[i],beta_covariates_sd_prior[i]);
  }

  for (i in 1:sum(D_cross[,1])) {
    cross_mu_fixed_coef[i] ~ normal(cross_mu_covariates_mean_prior[i],cross_mu_covariates_sd_prior[i]);
  }
  
  for (i in 1:sum(D_cross[,2])) {
    cross_rho_fixed_coef[i] ~ normal(cross_rho_covariates_mean_prior[i],cross_rho_covariates_sd_prior[i]);
  }

  if (prior_sim == 0) {
    for (k in 1:N_obs) {
      y_obs[k] ~ categorical_logit(x_beta[obs_idx[k]]');
    }
  }
  
}

generated quantities{
  int y_tilde[N];
  cov_matrix[2*T] Sigma;
  corr_matrix[2*T] Corr;
  matrix[N,K] x_beta;
  matrix[n,2*T] C;
  vector[T] PS_mu; //post sweep mu for identifiability
  vector[T] A_bar; //average of Ai
  vector[T] B_bar; //average of Bi

  C = (diag_pre_multiply(sigma, L_corr) * z)';
  for (t in 1:T){
    A_bar[t] = mean(C[,1 + 2 * (t - 1)]);
    B_bar[t] = mean(C[,2 + 2 *(t - 1)]);
    PS_mu[t] = mu[t] + A_bar[t] + B_bar[t];
  }
  Sigma = diag_pre_multiply(sigma, L_corr) * diag_pre_multiply(sigma, L_corr)';
  Corr = multiply_lower_tri_self_transpose(L_corr);

  // start of the x_beta calculation
  
  {
    
    // matrix[n,2*T] C; // for each actor, there are 2 * T number of random actor effects (two per network)
    matrix[n,T] alpha;
    matrix[n,T] beta;
    real mu_ij = 0;
    real mu_ji = 0;
    real rho_ij = 0;
    real cross_mu_ij = 0;
    real cross_rho_ij = 0;
    vector[n] alpha_fixed = rep_vector(0,n);
    vector[n] beta_fixed = rep_vector(0,n);
  
    // C = (diag_pre_multiply(sigma, L_corr) * z)';
    for (t in 1:T){
      int idx_a[2] = find_start_end(D_within[,3],t);
      int idx_b[2] = find_start_end(D_within[,4],t);


      if (D_within[t,3] > 0) {
        alpha_fixed = alpha_covariates[,idx_a[1]:idx_a[2]] * alpha_fixed_coef[idx_a[1]:idx_a[2]];
      }
      
      if (D_within[t,4] > 0) {
        beta_fixed = beta_covariates[,idx_b[1]:idx_b[2]] * beta_fixed_coef[idx_b[1]:idx_b[2]];
      }
      
      alpha[,t] = C[,1 + 2 * (t - 1)] + alpha_fixed;
      beta[,t] = C[,2 + 2 *(t - 1)] + beta_fixed;
    }
    
    {int counter = 1;
    matrix[T,2] M; 
      for (i in 1:n) {
        for (j in i:n) {
          if (i == j) {
            continue;
          }
          if (counter )
          for (k in 1:K){
  
            real within_terms = 0;
            real cross_terms = 0;
            
            for (t in 1:T) {
              real nt = ceil(k / 4^(t - 1));
              real score = fmod(nt,4);
              int idx_mu[2] = find_start_end(D_within[,1],t);
              int idx_rho[2] = find_start_end(D_within[,2],t);

              M[t,] = [(score == 2 || score == 0),(score == 3 || score == 0)];
              
              mu_ij = D_within[t,1] > 0 ? dot_product(mu_covariates[i,j][idx_mu[1]:idx_mu[2]], mu_fixed_coef[idx_mu[1]:idx_mu[2]]) : 0;
              mu_ji = D_within[t,1] > 0 ? dot_product(mu_covariates[j,i][idx_mu[1]:idx_mu[2]], mu_fixed_coef[idx_mu[1]:idx_mu[2]]) : 0;
              rho_ij = D_within[t,2] > 0 ? dot_product(rho_covariates[i,j][idx_rho[1]:idx_rho[2]], rho_fixed_coef[idx_rho[1]:idx_rho[2]]) : 0;
            

            //print("rho dot:", dot_product(rho_covariates[i,j][idx_rho[1]:idx_rho[2]], rho_fixed_coef[idx_rho[1]:idx_rho[2]]));
            
              within_terms += M[t,1]*(alpha[i,t] + beta[j,t] + mu[t] + mu_ij); 
              within_terms += M[t,2]*(alpha[j,t] + beta[i,t] + mu[t] + mu_ji); 
              within_terms += M[t,1]*M[t,2]*(rho[t] + rho_ij);
              //print("reciprocated: ", M[t,1]*M[t,2]); 

              
            }
            
            for (h in 1:H){
              int net_a = layer_pairs[h,1];
              int net_b = layer_pairs[h,2];
              
              int idx_mu[2] = find_start_end(D_cross[,1],h);
              int idx_rho[2] = find_start_end(D_cross[,2],h);


              cross_mu_ij = D_cross[h,1] > 0 ? dot_product(cross_mu_covariates[i,j][idx_mu[1]:idx_mu[2]],cross_mu_fixed_coef[idx_mu[1]:idx_mu[2]]): 0;
              cross_rho_ij = D_cross[h,2] > 0 ? dot_product(cross_rho_covariates[i,j][idx_rho[1]:idx_rho[2]], cross_rho_fixed_coef[idx_rho[1]:idx_rho[2]]) : 0;
                
                
              cross_terms += (M[net_a,1] * M[net_b,1] + M[net_a,2] * M[net_b,2]) * (cross_mu[h] + cross_mu_ij);
              cross_terms += (M[net_a,1] * M[net_b,2] + M[net_a,2] * M[net_b,1]) * (cross_rho[h] + cross_rho_ij);
            }
            
            x_beta[counter,k] = within_terms + cross_terms;
           
          }
          counter += 1;
        }
      }
    }
  }
  // end of the x_beta calculation
  
  for (k in 1:N) {
    y_tilde[k] = categorical_logit_rng(x_beta[k]');
  }

}

