
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
