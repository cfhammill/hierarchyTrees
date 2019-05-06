functions {
  real tree_lpf(int N, int R, int NNodes
               , int[] ranint_sizes, int[] node_number
               , matrix model_matrix, int[,] ranint_matrix
                , real tau, real tau_shape, real tau_rate
               , matrix L_omega, real lkj_shape
               , vector pivec, vector p_pi_conc
               , vector b_fix
               , vector[] err
               , matrix ranints, vector sigma_ranints, real ranint_shape
               , vector[] b
               , real sigma_model, real sig_model_shape
               , vector y){
    vector[N] y_est;
    real logProb;

    logProb = 0;
    logProb = logProb + gamma_lpdf(tau | tau_shape, tau_rate);
    logProb = logProb + lkj_corr_cholesky_lpdf(L_omega | lkj_shape);
    logProb = logProb + dirichlet_lpdf(pivec | p_pi_conc);
    logProb = logProb + normal_lpdf(b_fix | 0, 1);
    
    for(i in 1:NNodes){ // Does this assume a balanced tree...?
      logProb = logProb + normal_lpdf(err[i] | 0, 1);
    }
    
    for(r in 1:R) // This should probably be tabulated for each obs
      for(i in 1:(ranint_sizes[r]))
        logProb = logProb + normal_lpdf(ranints[i,r] | 0, sigma_ranints[r]);

    //formerly cauchy_lpdf(sigma_ranints | 0, ranint_shape);
    logProb = logProb + normal_lpdf(sigma_ranints | 0, ranint_shape);
    

    for(i in 1:N){
      y_est[i] = model_matrix[i,] * b_fix + model_matrix[i,] * b[node_number[i]]; //y = xb (1x2 * 2x1)

      for(r in 1:R)
        y_est[i] = y_est[i] + ranints[ranint_matrix[i,r], r];

      logProb = logProb + normal_lpdf(y[i] - y_est[i] | 0, sigma_model);
    }

    // formerly cauchy_lpdf(sigma_model | 0, sig_model_shape);
    logProb = logProb + normal_lpdf(sigma_model | 0, sig_model_shape);

    return logProb;
  }

  vector tree_rng(real tau_shape
                  , real tau_rate
                  , real lkj_shape
                  , vector p_pi_conc
                  , int NNodes, int P, int N
                  , int R, int[] ranint_sizes
                  , int[] node_parent, int[] node_number
                  , matrix model_matrix, int[,] ranint_matrix
                  , int Ranint_max
                  , real ranint_shape
                  , real sig_model_shape){
    
      real sigma_model;
      matrix[P,P] L_omega;
      real tau; 
      vector[P] pivec;
      vector[P] b_fix;
      
      vector[P] err[NNodes]; 
      
      matrix[Ranint_max,R] ranints;
      vector[R] sigma_ranints;
      matrix[P,P] L_sigma;
      vector[P] b[NNodes];
      vector[P] L_sigvec;

      vector[N] y_est;

      for(n in 1:N)
        y_est[n] = 0;

      for(p in 1:P)
        b_fix[p] = normal_rng(0,1);
      
      sigma_model = fabs(normal_rng(0, sig_model_shape)); //cauchy
      tau = fabs(gamma_rng(tau_shape, tau_rate));

      pivec = dirichlet_rng(p_pi_conc);

      L_omega = lkj_corr_cholesky_rng(P, lkj_shape);
      L_sigvec = tau * sqrt(P * pivec);
      L_sigma = diag_pre_multiply(L_sigvec, L_omega);
      
      if(node_parent[1] != 0)
        reject("Not at root! Something is amiss");

      for(i in 1:NNodes){
        int parent = node_parent[i];
        for(p in 1:P)
          err[i][p] = normal_rng(0,1);
        
        if(parent == 0){
          b[i] = L_sigma * err[i];
        } else {
          b[i] = b[node_parent[i]] + L_sigma * err[i];
        }
      }

      for(r in 1:R){
        sigma_ranints[r] = fabs(normal_rng(0, ranint_shape)); //cauchy
        for(i in 1:(ranint_sizes[r]))
          ranints[i,r] = normal_rng(0, sigma_ranints[r]);
      }

      for(i in 1:N){
        y_est[i] = model_matrix[i,] * b_fix + model_matrix[i,] * b[node_number[i]];

        for(r in 1:R){
          y_est[i] = y_est[i] + ranints[ranint_matrix[i,r], r];
        }

        y_est[i] = y_est[i] + normal_rng(0, sigma_model);
      }
      
      return y_est;
  }


  real tree_exp_lpf(int N, int R, int NNodes
               , int[] ranint_sizes, int[] node_number
               , matrix model_matrix, int[,] ranint_matrix
                , real tau, real tau_shape, real tau_rate
               , matrix L_omega, real lkj_shape
               , vector pivec, vector p_pi_conc
               , vector b_fix
               , vector[] err
               , matrix ranints, vector sigma_ranints, real ranint_shape
               , vector[] b
               , real sigma_model, real sig_model_shape
               , vector y){
    vector[N] y_est;
    real logProb;

    logProb = 0;
    logProb = logProb + gamma_lpdf(tau | tau_shape, tau_rate);
    logProb = logProb + lkj_corr_cholesky_lpdf(L_omega | lkj_shape);
    logProb = logProb + dirichlet_lpdf(pivec | p_pi_conc);
    logProb = logProb + normal_lpdf(b_fix | 0, 1);
    
    for(i in 1:NNodes){ // Does this assume a balanced tree...?
      logProb = logProb + normal_lpdf(err[i] | 0, 1);
    }
    
    for(r in 1:R) // This should probably be tabulated for each obs
      for(i in 1:(ranint_sizes[r]))
        logProb = logProb + normal_lpdf(ranints[i,r] | 0, sigma_ranints[r]);

    //formerly cauchy_lpdf(sigma_ranints | 0, ranint_shape);
    logProb = logProb + normal_lpdf(sigma_ranints | 0, ranint_shape);
    

    for(i in 1:N){
      y_est[i] = model_matrix[i,] * b_fix + model_matrix[i,] * b[node_number[i]]; //y = xb (1x2 * 2x1)

      for(r in 1:R)
        y_est[i] = y_est[i] + ranints[ranint_matrix[i,r], r];

      logProb = logProb + normal_lpdf(y[i] - y_est[i] | 0, sigma_model);
    }

    // formerly cauchy_lpdf(sigma_model | 0, sig_model_shape);
    logProb = logProb + exponential_lpdf(sigma_model | 1);

    return logProb;
  }
}
