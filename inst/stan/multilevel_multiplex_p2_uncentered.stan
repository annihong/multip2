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
      
      //given the number of actors in a group, index i and index j, find the index of the flattened dyadic covariates for group l
      int find_dyadic_covar_l_idx(int n_l, int i, int j) {
         int res = (i - 1) * (n_l - 1);
         int r = j < i ? j : j - 1;
         res += r;
         return res;
       }
      }
data {
  int<lower=0> prior_sim; //boolean for whether to simulate from priors (1) or fit the actual model (0)
  int<lower=0> network_sim; //boolean for whether to draw dyadic outcomes from the posterior, 1 yes and 0 no
  int<lower=0> L; //number of groups/repeated obs of the multiplex network 
  int<lower=0> n[L]; //number of actors for each of the L groups
  int<lower=0> N[L]; // number of total dyads for each of the L groups
  int<lower=0> T; // number of layers in the networks
  int<lower=0> H; // number of pairs of layers of networks
  int<lower=0> layer_pairs[H != 0 ? H : 1, 2]; //a list that stores all the network indices of the pairs of layers
  // if H > 0 (at least two networks), then H rows in the layer pairs else 1 row

  //outcomes (not gonna deal with missing dyads rn) 
  int<lower=0> K; //number of outcomes on a dyad := 2^2T
  int<lower=1, upper=K> y_obs[sum(N)]; // y_obs[sum_l(n_l)]
  
  //adding covariates to multiplex p2
  int<lower=0> D_within[T, 4]; // X[t, 1:4]  of covariates for mu_t, rho_t, alpha_t, beta_t
  int<lower=0> D_cross[H, 2]; // X[h, 1:2] of covar for cross_mu_h, cross_rho_h
  int<lower=0> N_covar[L]; // N_covar[l] = n_l(n_l - 1) for each group l
  
  vector[sum(D_within[,1])] mu_covariates[sum(N_covar)]; // for mu  
  vector[sum(D_within[,2])] rho_covariates[sum(N_covar)]; // for rho
  matrix[sum(n),sum(D_within[,3])] alpha_covariates; // for alpha
  matrix[sum(n),sum(D_within[,4])] beta_covariates; // for beta

  
  vector[H != 0 ? sum(D_cross[,1]) : 0] cross_mu_covariates[sum(N_covar)];
  vector[H != 0 ? sum(D_cross[,2]) : 0] cross_rho_covariates[sum(N_covar)];


// group level covariates
int<lower=0> p_group; // total number of group level covariates
int<lower=0> D_group_within[2, T]; // number of group level covariates for each layer
int<lower=0> D_group_cross[2, H]; // number of group level covariates for each pair of layers
matrix[L, p_group] group_covariates; // group level covariates
int<lower=0, upper=p_group>  mu_group_covariates_idx[sum(D_group_within[1])]; // indices of mu group level covariates
int<lower=0, upper=p_group>  rho_group_covariates_idx[sum(D_group_within[2])]; // indices of rho group level covariates
int<lower=0, upper=p_group>  cross_mu_group_covariates_idx[H != 0 ? sum(D_group_cross[1]) : 0]; // indices of cross mu group level covariates
int<lower=0, upper=p_group>  cross_rho_group_covariates_idx[H != 0 ? sum(D_group_cross[2]) : 0]; // indices of cross rho group level covariates



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
  vector[sum(D_group_within[1,])] mu_group_covariates_sd_prior;
  vector[sum(D_group_within[1,])] mu_group_covariates_mean_prior;

  vector[sum(D_group_within[2,])] rho_group_covariates_sd_prior;
  vector[sum(D_group_within[2,])] rho_group_covariates_mean_prior;

  vector[sum(D_within[,3])] alpha_covariates_sd_prior;
  vector[sum(D_within[,3])] alpha_covariates_mean_prior;

  vector[sum(D_within[,4])] beta_covariates_sd_prior;
  vector[sum(D_within[,4])] beta_covariates_mean_prior;

  vector[H != 0 ? sum(D_cross[,1]) : 0] cross_mu_covariates_sd_prior;
  vector[H != 0 ? sum(D_cross[,1]) : 0] cross_mu_covariates_mean_prior;

  vector[H != 0 ? sum(D_cross[,2]) : 0] cross_rho_covariates_sd_prior;
  vector[H != 0 ? sum(D_cross[,2]) : 0] cross_rho_covariates_mean_prior;

  vector[H != 0 ? sum(D_group_cross[1,]) : 0] cross_mu_group_covariates_sd_prior;
  vector[H != 0 ? sum(D_group_cross[1,]) : 0] cross_mu_group_covariates_mean_prior;

  vector[H != 0 ? sum(D_group_cross[2,]) : 0] cross_rho_group_covariates_sd_prior;
  vector[H != 0 ? sum(D_group_cross[2,]) : 0] cross_rho_group_covariates_mean_prior;
}

parameters {
  //the following are the overall intercepts pooling all L networks
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

  //group level fixed effects:
  vector[sum(D_group_within[1])] mu_group_coef;   
  vector[sum(D_group_within[2])] rho_group_coef; 
  vector[H != 0 ? sum(D_group_cross[1]) : 0] cross_mu_group_coef; 
  vector[H != 0 ? sum(D_group_cross[2]) : 0] cross_rho_group_coef;
  
  //random actor effects:
  //cov_matrix[2*T] Sigma[L]; // cov matrix used to draw the random actor effects
  vector<lower=0>[2*T] sigma;
  cholesky_factor_corr[2*T] L_corr;
  matrix[2*T, sum(n)] z; // z is flattend on n_l, so that different group can have different number of actors

  // random within layer effects (the 2*T will change when we include more within layer random effects)
  //cov_matrix[2*L] Sigma_within[T]; 
  vector<lower=0>[2] sigma_within[T];
  cholesky_factor_corr[2] L_corr_within[T]; 
  vector[L] z_mu[T];
  vector[L] z_rho[T];
  //z_within[T] = [z_mu, z_rho]

  // random cross layer effects (the 2*T will change when we include more cross layer random effects)
  //cov_matrix[2*L] Sigma_cross[H]; 
  vector<lower=0>[2] sigma_cross[H];
  cholesky_factor_corr[2] L_corr_cross[H];
  vector[L] z_cross_mu[H];
  vector[L] z_cross_rho[H];
  //z_cross[H] = [z_cross_mu, z_cross_rho]





}

transformed parameters{
    matrix[sum(N),K] x_beta;
    matrix[L, T] mu_random;
    matrix[L, T] rho_random;
    matrix[L, H] cross_mu_random;
    matrix[L, H] cross_rho_random;

    matrix[L, T] mu_group_fixed = rep_matrix(0, L, T);
    matrix[L, T] rho_group_fixed = rep_matrix(0, L, T);
    matrix[L, H] cross_mu_group_fixed = rep_matrix(0, L, H);
    matrix[L, H] cross_rho_group_fixed = rep_matrix(0, L, H);


  {

    for (t in 1:T) {
      matrix[L,2] C_within;
      matrix[2,L] z_within_t = append_col(z_mu[t], z_rho[t])';
      C_within = (diag_pre_multiply(sigma_within[t], L_corr_within[t]) * z_within_t)';
      mu_random[,t] = C_within[,1] + mu[t];
      rho_random[,t] = C_within[,2] + rho[t];
      // print("mu_group_covariates_idx: ", mu_group_covariates_idx);
      // print("rho_group_covariates_idx: ", rho_group_covariates_idx);      
      
      if (D_group_within[1,t] > 0) {
        int idx_mu_group[2] = find_start_end(D_group_within[1],t);
        mu_group_fixed[,t] = (group_covariates[,mu_group_covariates_idx[idx_mu_group[1]:idx_mu_group[2]]] * mu_group_coef[idx_mu_group[1]:idx_mu_group[2]]);
      }
        
      if (D_group_within[2,t] > 0) {
        int idx_rho_group[2] = find_start_end(D_group_within[2],t);
        rho_group_fixed[,t] = (group_covariates[,rho_group_covariates_idx[idx_rho_group[1]:idx_rho_group[2]]] * rho_group_coef[idx_rho_group[1]:idx_rho_group[2]]);
        }  
            // print("mu_group_fixed: ", mu_group_fixed);
      // print("rho_group_fixed: ", rho_group_fixed);
      //print(C_within);
    }

    // print("mu_random: ", mu_random);
    // print("rho_random: ", rho_random);

    for (h in 1:H) {
      matrix[L,2] C_cross;
      matrix[2,L] z_cross_h = append_col(z_cross_mu[h], z_cross_rho[h])';
      C_cross = (diag_pre_multiply(sigma_cross[h], L_corr_cross[h]) * z_cross_h)';
      cross_mu_random[,h] = C_cross[,1] + cross_mu[h];
      cross_rho_random[,h] = C_cross[,2] + cross_rho[h];

      if (D_group_cross[1,h] > 0) {
        int idx_c_mu_group[2] = find_start_end(D_group_cross[1],h);
        cross_mu_group_fixed[,h] = (group_covariates[,cross_mu_group_covariates_idx[idx_c_mu_group[1]:idx_c_mu_group[2]]] * cross_mu_group_coef[idx_c_mu_group[1]:idx_c_mu_group[2]]);
      }

      if (D_group_cross[2,h] > 0) {
        int idx_c_rho_group[2] = find_start_end(D_group_cross[2],h);
        cross_rho_group_fixed[,h] = (group_covariates[,cross_rho_group_covariates_idx[idx_c_rho_group[1]:idx_c_rho_group[2]]] * cross_rho_group_coef[idx_c_rho_group[1]:idx_c_rho_group[2]]);
      }
    }

  }

  {   // start of the x_beta calculation
  // start of looping over L groups of T-plex networks 
    for (l in 1:L) {
      int n_l = n[l];
      matrix[n[l],2*T] C; // for each actor, there are 2 * T number of random actor effects (two per network)
      matrix[n[l],T] alpha;
      matrix[n[l],T] beta;
      real mu_ij = 0;
      real mu_ji = 0;
      real rho_ij = 0;
      real cross_mu_ij = 0;
      real cross_rho_ij = 0;
      vector[n[l]] alpha_fixed = rep_vector(0,n[l]);
      vector[n[l]] beta_fixed = rep_vector(0,n[l]);
      int idx_nl[2] = find_start_end(n,l);
      int idx_N_covar_l[2] = find_start_end(N_covar,l);
      C = (diag_pre_multiply(sigma, L_corr) * z[,idx_nl[1]:idx_nl[2]])';
      

      // start of looping over T layers to update random actor effects:
      for (t in 1:T){
        int idx_a[2] = find_start_end(D_within[,3],t);
        int idx_b[2] = find_start_end(D_within[,4],t);


        if (D_within[t,3] > 0) {
          alpha_fixed = alpha_covariates[idx_nl[1]:idx_nl[2],idx_a[1]:idx_a[2]] * alpha_fixed_coef[idx_a[1]:idx_a[2]];
        }
        
        if (D_within[t,4] > 0) {
          beta_fixed = beta_covariates[idx_nl[1]:idx_nl[2] ,idx_b[1]:idx_b[2]] * beta_fixed_coef[idx_b[1]:idx_b[2]];
        }
        
        alpha[,t] = C[,1 + 2 * (t - 1)] + alpha_fixed;
        beta[,t] = C[,2 + 2 *(t - 1)] + beta_fixed;
      }
      // end of looping over T layers to update random actor effects:
      
      {

      int counter = sum(N[1:(l-1)]) + 1; //the the counter to start on the first dyad of the group
      // print("counter: ", counter);
      // print("l: ", l);

      matrix[T,2] M; 

        for (i in 1:n_l) {
          for (j in i:n_l) {
            if (i == j) {
              continue;
            }
            //if (counter ) ??

            // start of looping through all K outcomes:
            for (k in 1:K){
    
              real within_terms = 0;
              real cross_terms = 0;
              int idx_ij_dyad = idx_N_covar_l[1] + find_dyadic_covar_l_idx(n_l, i, j) - 1;
              int idx_ji_dyad = idx_N_covar_l[1] + find_dyadic_covar_l_idx(n_l, j, i) - 1;
              
              for (t in 1:T) {
                real nt = ceil(k / 4^(t - 1));
                real score = fmod(nt,4);


                M[t,] = [(score == 2 || score == 0),(score == 3 || score == 0)];
                
                if (D_within[t,1] > 0) {
                    int idx_mu[2] = find_start_end(D_within[,1],t);
                    mu_ij = dot_product(mu_covariates[idx_ij_dyad][idx_mu[1]:idx_mu[2]], mu_fixed_coef[idx_mu[1]:idx_mu[2]]);
                    mu_ji = dot_product(mu_covariates[idx_ji_dyad][idx_mu[1]:idx_mu[2]], mu_fixed_coef[idx_mu[1]:idx_mu[2]]);
                }

                if (D_within[t,2] > 0) {
                  int idx_rho[2] = find_start_end(D_within[,2],t);
                  rho_ij = dot_product(rho_covariates[idx_ij_dyad][idx_rho[1]:idx_rho[2]], rho_fixed_coef[idx_rho[1]:idx_rho[2]]);
                }
              //print("rho dot:", dot_product(rho_covariates[i,j][idx_rho[1]:idx_rho[2]], rho_fixed_coef[idx_rho[1]:idx_rho[2]]));
              
                within_terms += M[t,1]*(alpha[i,t] + beta[j,t] + mu_random[l, t] + mu_group_fixed[l,t] +  mu_ij); 
                within_terms += M[t,2]*(alpha[j,t] + beta[i,t] + mu_random[l, t] + mu_group_fixed[l,t] + mu_ji); 
                within_terms += M[t,1]*M[t,2]*(rho_random[l,t] + rho_group_fixed[l,t] + rho_ij);
                //print("reciprocated: ", M[t,1]*M[t,2]); 

                
              }
              
              for (h in 1:H){
                int net_a = layer_pairs[h,1];
                int net_b = layer_pairs[h,2];
                
                if (D_cross[h,1] > 0) {
                  int idx_c_mu[2] = find_start_end(D_cross[,1],h);
                  cross_mu_ij = dot_product(cross_mu_covariates[idx_ij_dyad][idx_c_mu[1]:idx_c_mu[2]],cross_mu_fixed_coef[idx_c_mu[1]:idx_c_mu[2]]);
                }

                if (D_cross[h,2] > 0) {
                  int idx_c_rho[2] = find_start_end(D_cross[,2],h);
                  cross_rho_ij = dot_product(cross_rho_covariates[idx_ij_dyad][idx_c_rho[1]:idx_c_rho[2]], cross_rho_fixed_coef[idx_c_rho[1]:idx_c_rho[2]]);
                }
                  
                cross_terms += (M[net_a,1] * M[net_b,1] + M[net_a,2] * M[net_b,2]) * (cross_mu_random[l, h] + cross_mu_group_fixed[l,h] + cross_mu_ij);
                cross_terms += (M[net_a,1] * M[net_b,2] + M[net_a,2] * M[net_b,1]) * (cross_rho_random[l, h]+ cross_rho_group_fixed[l,h] + cross_rho_ij);
              }
              
              x_beta[counter,k] = within_terms + cross_terms;
            
            }
            // end of looping through all K outcomes:
            counter += 1;
          }
        }
      }
    }
      // end of looping over L groups of T-plex networks 
  }
  

// print("cross_mu_random: ", cross_mu_random);
// print("cross_rho_random: ", cross_rho_random);
}


model {
  // end of the x_beta calculation


  
  for (t in 1:T) {
    mu[t] ~ normal(mu_mean_prior[t],mu_sd_prior[t]);
    rho[t] ~ normal(rho_mean_prior[t],rho_sd_prior[t]);
    z_mu[t] ~ normal(0, 1);
    z_rho[t] ~ normal(0, 1);
    L_corr_within[t] ~ lkj_corr_cholesky(LJK_eta_prior);
    sigma_within[t] ~ inv_gamma(scale_alpha_prior,scale_beta_prior);
  }



  for (h in 1:H) {
    cross_mu[h] ~ normal(cross_mu_mean_prior[h],cross_mu_sd_prior[h]);
    cross_rho[h] ~ normal(cross_rho_mean_prior[h],cross_rho_sd_prior[h]);
    z_cross_mu[h] ~ normal(0, 1);
    z_cross_rho[h] ~ normal(0, 1);
    L_corr_cross[h] ~ lkj_corr_cholesky(LJK_eta_prior);
    sigma_cross[h] ~ inv_gamma(scale_alpha_prior,scale_beta_prior);
  }

    //for random actor effects:
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

  //priors for the group level coefficients

  for (i in 1:sum(D_group_within[1])) {
    mu_group_coef[i] ~ normal(mu_group_covariates_mean_prior[i],mu_group_covariates_sd_prior[i]);
  }
  for (i in 1:sum(D_group_within[2])) {
    rho_group_coef[i] ~ normal(rho_group_covariates_mean_prior[i],rho_group_covariates_sd_prior[i]);
  }
  for (i in 1:sum(D_group_cross[1])) {
    cross_mu_group_coef[i] ~ normal(cross_mu_group_covariates_mean_prior[i],cross_mu_group_covariates_sd_prior[i]);
  }
  for (i in 1:sum(D_group_cross[2])) {
    cross_rho_group_coef[i] ~ normal(cross_rho_group_covariates_mean_prior[i],cross_rho_group_covariates_sd_prior[i]);
  }


  if (prior_sim == 0) {
    for (dyad_idx in 1:sum(N)) {
      y_obs[dyad_idx] ~ categorical_logit(x_beta[dyad_idx]');
    }
  }
  
}

generated quantities{
  int y_tilde[sum(N)];
  cov_matrix[2*T] Sigma;
  corr_matrix[2*T] Corr;
  cov_matrix[2] Sigma_within[T];
  corr_matrix[2] Corr_within[T];
  cov_matrix[2] Sigma_cross[H];
  corr_matrix[2] Corr_cross[H];
  vector[T] mu_PS; // vector of size T for within-network density
  vector[T] rho_PS; // vector of size T for within-network reciprocity
  vector[H] cross_mu_PS; // size H storing cross-network density for each pair of network
  vector[H] cross_rho_PS; // size H storing cross-network reciprocity for each pair of network
  matrix[L, T] mu_random_PS;
  matrix[L, T] rho_random_PS;
  matrix[L, H] cross_mu_random_PS;
  matrix[L, H] cross_rho_random_PS;
  real A_bar[T];
  real B_bar[T];
  real mu_random_bar[T];
  real rho_random_bar[T];
  real cross_mu_random_bar[H];
  real cross_rho_random_bar[H];

  Sigma = diag_pre_multiply(sigma, L_corr) * diag_pre_multiply(sigma, L_corr)';
  Corr = multiply_lower_tri_self_transpose(L_corr);




  {
    matrix[L,T] A_bar_temp;
    matrix[L,T] B_bar_temp;
    for (l in 1:L) {
      matrix[n[l],2*T] C; // for each actor, there are 2 * T number of random actor effects (two per network)
      int idx_nl[2] = find_start_end(n,l);
      C = (diag_pre_multiply(sigma, L_corr) * z[,idx_nl[1]:idx_nl[2]])';
      for (t in 1:T){
        A_bar_temp[l,t] = mean(C[,1 + 2 * (t - 1)]);
        B_bar_temp[l,t] = mean(C[,2 + 2 *(t - 1)]);
      }
    }
  
    
    for (t in 1:T) {
      // matrix[2,L] C_within;
      // matrix[2,L] z_within_t = append_col(z_mu[t], z_rho[t])';

      Sigma_within[t] = diag_pre_multiply(sigma_within[t], L_corr_within[t]) * diag_pre_multiply(sigma_within[t], L_corr_within[t])';
      Corr_within[t] = multiply_lower_tri_self_transpose(L_corr_within[t]);

      // C_within = (diag_pre_multiply(sigma_within[t], L_corr_within[t]) * z_within_t);

      mu_random_bar[t] = mean(mu_random[,t] - mu[t]);
      rho_random_bar[t] = mean(rho_random[,t] - rho[t]);
      A_bar[t] = mean(A_bar_temp[,t]);
      B_bar[t] = mean(B_bar_temp[,t]);
      mu_PS[t] = mu[t] + mu_random_bar[t] + A_bar[t] + B_bar[t];
      rho_PS[t] = rho[t] + rho_random_bar[t];
      mu_random_PS[,t] = mu_random[,t] - mu_random_bar[t] - A_bar[t] - B_bar[t] + A_bar_temp[,t] + B_bar_temp[,t];
      rho_random_PS[,t] = rho_random[,t] - rho_random_bar[t];
      
    }

    for (h in 1:H) {
      // matrix[2,L] C_cross;
      // matrix[2,L] z_cross_h = append_col(z_cross_mu[h], z_cross_rho[h])';

      Sigma_cross[h] = diag_pre_multiply(sigma_cross[h], L_corr_cross[h]) * diag_pre_multiply(sigma_cross[h], L_corr_cross[h])';
      Corr_cross[h] = multiply_lower_tri_self_transpose(L_corr_cross[h]);
      // C_cross = (diag_pre_multiply(sigma_cross[h], L_corr_cross[h]) * z_cross_h);
      cross_mu_random_bar[h] = mean(cross_mu_random[,h] - cross_mu[h]);
      cross_rho_random_bar[h] = mean(cross_rho_random[,h] - cross_rho[h]);

      cross_mu_PS[h] = cross_mu[h] + cross_mu_random_bar[h];
      cross_rho_PS[h] = cross_rho[h] + cross_rho_random_bar[h];
      cross_mu_random_PS[,h] = cross_mu_random[,h] - cross_mu_random_bar[h];s
      cross_rho_random_PS[,h] = cross_rho_random[,h] - cross_rho_random_bar[h];
    }
  }



  // start of the x_beta calculation
  
  if (network_sim == 1) { //only draw from the posterior when it is yes
    for (k in 1:sum(N)) {
      y_tilde[k] = categorical_logit_rng(x_beta[k]');
    }
  }

}

