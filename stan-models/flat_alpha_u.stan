data {
  int <lower=1> M; // number of countries
  int <lower=1> P; // number of covariates
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N[M]; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  int cases[N2,M]; // reported cases
  int deaths[N2, M]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  matrix[N2, M] f; // h * s
  matrix[N2, P] X[M]; // features matrix
  int EpidemicStart[M];
  real pop[M];
  real SI[N2]; // fixed pre-calculated SI using emprical data from Neil
}

transformed data {
  vector[N2] SI_rev; // SI in reverse order
  vector[N2] f_rev[M]; // f in reversed order
  matrix[N2, M] zero_matrix = rep_matrix(0, N2, M);
  
  for(i in 1:N2)
    SI_rev[i] = SI[N2-i+1];
    
  for(m in 1:M){
    for(i in 1:N2) {
     f_rev[m, i] = f[N2-i+1,m];
    }
  }
}


parameters {
  real<lower=0> mu[M]; // intercept for Rt
  vector<lower=0, upper = 1.2>[P] alpha;
  real<lower=0> gamma;
  vector[M] lockdown;
  real<lower=0> kappa;
  real<lower=0> y[M];
  real<lower=0> phi;
  real<lower=0> tau_unit;
  vector<lower=0>[M] ifr_noise;
  real log_infecteds_multiplier; // We will _divide_ ifr by exp(log_infecteds_multiplier)
  real <lower = 0, upper = 1> u; // Proportion of population that is isolated after lockdown
}

transformed parameters {
    matrix[N2, M] prediction = zero_matrix;
    matrix[N2, M] E_deaths  = zero_matrix;
    matrix[N2, M] Rt = zero_matrix;
    matrix[N2, M] Rt_adj = Rt;
    // The IFR calculated by Imperial group will be multiplied by the random ifr_factor
    vector<lower = 0>[M] ifr_factor = ifr_noise / exp(log_infecteds_multiplier); 
    // Number of initial infecteds should also be multiplied by infectes_multiplier,
    // hence tau ~ exponential(0.03) * infecteds_multiplier
    real<lower = 0> tau = tau_unit / 0.03 * exp(log_infecteds_multiplier); 
    
    {
      matrix[N2,M] cumm_sum = zero_matrix;
      matrix[N2,M] cumm_sum_lockdown = zero_matrix;
      for (m in 1:M){
        prediction[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
        cumm_sum[2:N0,m] = cumulative_sum(prediction[2:N0,m]);
        cumm_sum_lockdown[2:N0, m] = cumm_sum[2:N0, m];
        
        Rt[,m] = mu[m] * exp(-X[m] * alpha - X[m][,5] * lockdown[m]);
        Rt_adj[1:N0,m] = Rt[1:N0,m];
        for (i in (N0+1):N2) {
          real convolution = dot_product(sub_col(prediction, 1, m, i-1), tail(SI_rev, i-1));
          cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i-1,m];
          cumm_sum_lockdown[i,m] = cumm_sum[i-1,m] + prediction[i-1,m] * (1 - X[m][i, 5]);
          Rt_adj[i,m] = ((pop[m] - u * (pop[m] - cumm_sum_lockdown[i, m]) * X[m][i, 5] - cumm_sum[i,m]) / pop[m]) * Rt[i,m];
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
  for (m in 1:M){
      y[m] ~ exponential(1/tau);
  }
  gamma ~ normal(0,.2);
  lockdown ~ normal(0,gamma);
  phi ~ normal(0,5);
  kappa ~ normal(0,0.5);
  mu ~ normal(3.28, kappa); // citation: https://academic.oup.com/jtm/article/27/2/taaa021/5735319
  ifr_noise ~ normal(1,0.1);
  log_infecteds_multiplier ~ normal(0, 2);
  for(m in 1:M){
    deaths[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:N[m], m], phi);
   }
}

generated quantities {
    matrix[N2, M] prediction0 = zero_matrix;
    matrix[N2, M] E_deaths0  = zero_matrix;
    
    {
      matrix[N2,M] cumm_sum0 = zero_matrix;
      for (m in 1:M){
         for (i in 2:N0){
          cumm_sum0[i,m] = cumm_sum0[i-1,m] + y[m]; 
        }
        prediction0[1:N0,m] = rep_vector(y[m],N0); 
        for (i in (N0+1):N2) {
          real convolution0 = dot_product(sub_col(prediction0, 1, m, i-1), tail(SI_rev, i-1));
          cumm_sum0[i,m] = cumm_sum0[i-1,m] + prediction0[i-1,m];
          prediction0[i, m] = ((pop[m]-cumm_sum0[i,m]) / pop[m]) * mu[m] * convolution0;
        }
        E_deaths0[1, m]= 1e-15 * prediction0[1,m];
        for (i in 2:N2){
          E_deaths0[i,m] = ifr_factor[m] * dot_product(sub_col(prediction0, 1, m, i-1), tail(f_rev[m], i-1));
        }
      }
    }
}

