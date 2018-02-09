data{
  int N;
  int NNodes;
  int NParents;
  int node_parent[NNodes]; // Node to parent number mapping
  int P;
  int R;
  int Ranint_max;
  int ranint_sizes[R];
  matrix[N,P] model_matrix;
  int ranint_matrix[N,R];
  int node_number[N];
  vector[N] y;
  real lkj_shape;
  real root_shape;
  real sig_shape;
  real ranint_shape;
  real tau_shape;
  real pi_conc;
  real sig_model_shape;
}

transformed data{
  vector[P] p_zeros;
  vector[P] p_pi_conc;
  
  for(i in 1:P){
    p_zeros[i] = 0;
    p_pi_conc[i] = pi_conc;
  }
}

parameters{
  real<lower=0> sigma_model;

  cholesky_factor_corr[P] L_omega;
  real<lower=0> tau; 
  simplex[P] pivec; 

  vector[P] err[NNodes]; 
  
  matrix[Ranint_max,R] ranints;
  vector<lower=0>[R] sigma_ranints;
}
         
transformed parameters{
  cholesky_factor_cov[P] L_sigma;
  vector[P] b[NNodes];
  vector<lower=0>[P] L_sigvec;

  
  L_sigvec = tau * sqrt(P * pivec);
  L_sigma = diag_pre_multiply(L_sigvec, L_omega);

  if(node_parent[1] != 0)
    reject("Not at root! Something is amiss");

  for(i in 1:NNodes){
    int parent = node_parent[i];
    if(parent == 0){
      b[i] = L_sigma * err[i];
    } else {
      b[i] = b[node_parent[i]] + L_sigma * err[i];
    }
  }
}
  
         
model{
  vector[N] y_est;
  target += gamma_lpdf(tau | tau_shape, 1);
  target += lkj_corr_cholesky_lpdf(L_omega | lkj_shape);
  target += dirichlet_lpdf(pivec | p_pi_conc);     

  for(i in 1:NNodes){ // Does this assume a balanced tree...?
    target += normal_lpdf(err[i] | 0, 1);
  }
  
  for(r in 1:R) // This should probably be tabulated for each obs
    for(i in 1:(ranint_sizes[r]))
      target += normal_lpdf(ranints[i,r] | 0, sigma_ranints[r]);

  target += cauchy_lpdf(sigma_ranints | 0, ranint_shape);

  for(i in 1:N){
    y_est[i] = model_matrix[i,] * b[node_number[i]]; //y = xb (1x2 * 2x1)

    for(r in 1:R)
      y_est[i] = y_est[i] + ranints[ranint_matrix[i,r], r];

    target += normal_lpdf(y[i] - y_est[i] | 0, sigma_model);
  }
  
  target += cauchy_lpdf(sigma_model | 0, sig_model_shape);
}

generated quantities {
  vector[N] logLik; 

  {
    vector[N] y_pred;
    for(i in 1:N){
   
      y_pred[i] = model_matrix[i,] * b[node_number[i]]; //y = xb (1x2 * 2x1)

      for(r in 1:R)
	y_pred[i] = y_pred[i] + ranints[ranint_matrix[i,r], r];
    
      logLik[i] = normal_lpdf(y[i] | y_pred[i], sigma_model);
    }
  }
  
}
