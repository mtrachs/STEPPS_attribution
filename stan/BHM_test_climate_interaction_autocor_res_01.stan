
data {

int S; // number of grid cells
int time_slice; // number of time points
int K; // number of taxa

matrix<lower=0,upper=1>[K,time_slice] y[S]; // vegetation reconstruction, currently has format different from format defined below!!!!
vector[time_slice] M; //moisture
vector[time_slice] T; //temperature
vector[time_slice] M2; //moisture
vector[time_slice] T2; //temperature
vector[time_slice] TM;
matrix[S,S] d;
vector[S] coord_x;
vector[S] coord_y;
vector[time_slice] TM_01; //interaction temperature 

}

transformed data {
  
  vector[S] zeros;
  vector[K] zeros_k;
  vector[K] var_k;
  //vector[time_slice] T_sq; 
  //vector[time_slice] M_sq;
  
  
  
  for(i in 1:S) zeros[i] = 0;
  for(k in 1:K) zeros_k[k] = 0;
  for(k in 1:K) var_k[k] = 10;
  //for(t in 1:time_slice) T_sq[t] = T[t]*T[t];
  //for(t in 1:time_slice) M_sq[t] = M[t]*M[t];

  
  //probably add spatial data here 

}

parameters{
  
  //parameters of dirichlet distribution and its components
    vector<lower=-10,upper=10>[K] alpha_zero_k;
//alpha_T_K_t; // interaction term 
   /* matrix<lower=-10,upper=10>[time_slice,K] alpha_one;
    matrix<lower=-10,upper=10>[time_slice,K] alpha_two;
    matrix<lower=-10,upper=10>[time_slice,K] alpha_three;
    matrix<lower=-10,upper=10>[time_slice,K] alpha_four;
*/
    #vector<lower=-10,upper=10>[K] alpha_one;
    #vector<lower=-10,upper=10>[K] alpha_two;
    vector<lower=-10,upper=10>[K] alpha_three;
    vector<lower=-10,upper=10>[K] alpha_four;
    vector<lower=-10,upper=10>[K] alpha_seven;
    //vector<lower=-10,upper=10>[K] alpha_six;

  /*vector[K] alpha_zero_k;
  matrix[time_slice,K] alpha_M_t;
  matrix[time_slice,K] alpha_T_t;
  //alpha_T_K_t; // interaction term 
  matrix[time_slice,K] alpha_one;
  matrix[time_slice,K] alpha_two;
*/  //vector alpha_three[K];
  
  //matrix<lower=-10,upper=10>[time_slice,K] eps_t_k;
  //matrix<lower=-10,upper=10>[time_slice,K] eps_s_t_k[S];
  //cov_matrix[S] eps_s_k[K];// should Stan estimate this or do we define it?
 vector<lower=-10,upper=10>[S] mu_s_k[K];
  vector<lower=1e-4,upper=10>[K] eta;
  vector<lower=1e-4,upper=10>[K] rho;
  
}

/*transformed parameters {

as log to ensure alpha is positive 
}*/


model {
  vector[K] y_comp;//matrix[S,time_slice] alpha_s_t_k[K];
  matrix[time_slice,K] clim_rel;
  vector[K] alpha_s_t_k;
  matrix[S,S] C_s[K];
  matrix[S,S] C_s_L[K];
  
  
  alpha_zero_k ~ normal(0,10);
  #alpha_one ~ multi_normal(zeros_k,diag_matrix(var_k));
  #alpha_two ~ multi_normal(zeros_k,diag_matrix(var_k));
  alpha_three ~ multi_normal(zeros_k,diag_matrix(var_k));
  alpha_four ~ multi_normal(zeros_k,diag_matrix(var_k));
  alpha_seven ~ multi_normal(zeros_k,diag_matrix(var_k));
  //alpha_six ~ multi_normal(zeros_k,diag_matrix(var_k));
  
  
    
    for(k in 1:K){
      C_s[k]   = exp(-d/rho[k]);
      C_s_L[k] = cholesky_decompose(C_s[k]);
      mu_s_k[k] ~ multi_normal_cholesky(zeros,eta[k]*C_s_L[k]);
     
    }
  

  
    for(s in 1:S){
      for(t in 1:time_slice){
        for(k in 1:K){
  // perhaps omit eps_t_k[t,k] in a first pass
          
        clim_rel[t,k] = alpha_zero_k[k]  + alpha_three[k]*T[t] +  alpha_four[k]*M[t] + alpha_seven[k]*TM_01[t];                  alpha_s_t_k[k] = clim_rel[t,k] + mu_s_k[k,s];  
        y_comp[k] = y[s,k,t];
        }
      // specifiy the likelihood
        //print(alpha_s_t_k);
        y_comp ~ dirichlet(exp(alpha_s_t_k));
      }
    }
}


generated quantities {
  matrix[S,time_slice] log_lik;
  {
  vector[K] y_comp;//matrix[S,time_slice] alpha_s_t_k[K];
  matrix[time_slice,K] clim_rel;
  vector[K] alpha_s_t_k;

   for(s in 1:S){
      for(t in 1:time_slice){
        for(k in 1:K){
  // perhaps omit eps_t_k[t,k] in a first pass
          
          clim_rel[t,k] = alpha_zero_k[k]  + alpha_three[k]*T[t] +  alpha_four[k]*M[t] + alpha_seven[k]*TM_01[t];                alpha_s_t_k[k] = clim_rel[t,k] + mu_s_k[k,s];  
          y_comp[k] = y[s,k,t];
        }
       // specifiy the likelihood
        //print(alpha_s_t_k);
       log_lik[s,t] = dirichlet_lpdf(y_comp|exp(alpha_s_t_k));
      }
    }
} 
}