#include /diffusion_tree_functions.stan

data{
  // Data dimensions
  int N;                    // Number of observations
  int P;                    // Number of node dependent effects
  int R;                    // Number of random intercepts
  int NNodes;               // Number of nodes
  int ranint_sizes[R];      // Number of levels for each random intercept
  int Ranint_max;           // Max levels for the random intercepts
  
  // Observation data
  int node_number[N];       // Which node each volume is
  int node_parent[NNodes];  // Node to parent number mapping
  matrix[N,P] model_matrix; // Model matrix for node dependent effects
  int ranint_matrix[N,R];   // The random intercept levels for each obs

  // Prior shapes
  real ranint_shape;        // Cauchy shape for random effect variance
  real lkj_shape;           // Shape for the LKJ for within node eff cov          
  real tau_shape;           // Cauchy shape for LKJ tau
  real tau_rate;
  real pi_conc;             // Concentration param for LKJ
  real sig_model_shape;     // Cauchy shape for the model variance
}

transformed data{
  vector[P] p_zeros;
  vector[P] p_pi_conc;
  
  for(i in 1:P){
    p_zeros[i] = 0;
    p_pi_conc[i] = pi_conc;
  }
}

generated quantities {
  // I feel bad including a whole chunk of program here, but it is too unweildy
  // to return a massive vector of results and then unpack it from the rng
  // function.

#include /diffusion_tree_rng.stan
}
