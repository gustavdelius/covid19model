data {
  int <lower=1> M; // number of countries
  int <lower=1> P; // number of covariates for full pooling (global effects)
  int <lower=1> P_partial; // number of covariates for partial pooling (state-level effects)
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N[M]; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  int deaths[N2, M]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  matrix[N2, M] f; // ifr
  matrix[N2, P] X[M];
  matrix[N2, P_partial] X_partial[M];
  int EpidemicStart[M];
  real pop[M];
  real SI[N2]; // fixed SI using empirical data
}

transformed data {
  vector[N2] SI_rev; // SI in reverse order
  vector[N2] f_rev[M]; // f in reversed order
  
  for(i in 1:N2)
    SI_rev[i] = SI[N2-i+1];
    
  for(m in 1:M){
    for(i in 1:N2) {
     f_rev[m, i] = f[N2-i+1,m];
    }
  }
}


parameters {
  real<lower=0> R0[M]; 
  vector[P] alpha; 
  vector[P_partial] alpha_state[M];
  real<lower=0> gamma;
  real<lower=0> kappa;
  real<lower=0> y[M];
  real<lower=0> phi;
  real<lower=0> tau_unit;
  vector<lower=0>[M] ifr_noise;
  real mu;
}

transformed parameters {
    matrix[N2, M] prediction = rep_matrix(0,N2,M);
    matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
    matrix[N2, M] Rt = rep_matrix(0,N2,M);
    matrix[N2, M] Rt_adj = Rt;
    // The IFR calculated by Imperial group will be multiplied by the random ifr_factor
    vector<lower = 0>[M] ifr_factor = ifr_noise * exp(mu); 
    // Number of initial infecteds should also be multiplied by infectes_multiplier,
    // hence tau ~ exponential(0.03) * infecteds_multiplier
    real<lower = 0> tau = tau_unit / 0.03 / exp(mu); 
    
    {
      matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
      for (m in 1:M){
        prediction[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
        cumm_sum[2:N0,m] = cumulative_sum(prediction[2:N0,m]);
        
        Rt[,m] = R0[m] * 2 * inv_logit(-X[m] * alpha - X_partial[m] * alpha_state[m]);
        Rt_adj[1:N0,m] = Rt[1:N0,m];
        
        for (i in (N0+1):N2) {
          real convolution = dot_product(sub_col(prediction, 1, m, i-1), tail(SI_rev, i-1));
          
          cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i-1,m];
          Rt_adj[i,m] = ((pop[m]-cumm_sum[i,m]) / pop[m]) * Rt[i,m];
          prediction[i, m] = Rt_adj[i,m] * convolution;
        }
        
        E_deaths[1, m]= 1e-15 * prediction[1,m];
        for (i in 2:N2){
          E_deaths[i,m] = ifr_factor[m] * dot_product(sub_col(prediction, 1, m, i-1), tail(f_rev[m], i-1));
        }
      }
    }
}
model {
  tau_unit ~ exponential(1);
  gamma ~ normal(0,.5);
  for (m in 1:M) {
      y[m] ~ exponential(1/tau);
      alpha_state[m] ~ normal(0,gamma);
  }
  phi ~ normal(0,5);
  kappa ~ normal(0,0.5);
  R0 ~ normal(3.28, kappa); // citation: https://academic.oup.com/jtm/article/27/2/taaa021/5735319
  alpha ~ normal(0,0.5);
  ifr_noise ~ normal(1,0.1);
  mu ~ normal(0, 2);

  for(m in 1:M){
    deaths[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:N[m], m], phi);
   }
}
