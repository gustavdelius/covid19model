data {
  int<lower = 1> M; // number of countries
  int<lower = 1> N0; // number of days for which to impute infections
  int<lower = 1> N[M]; // number of days of observed data for country m. each entry must be <= N2
  int<lower = 1> N2; // number of days of observed data + # of days to forecast
  int cases[N2, M]; // number of reported cases -- the rows with i > N contain -1 and should be ignored
  int deaths[N2, M]; // number of reported deaths -- set to -1 for i > N
  matrix[N2, M] f; // mean daily mortality rate of infecteds
  matrix[N2, M] covariate1;
  matrix[N2, M] covariate2;
  matrix[N2, M] covariate3;
  matrix[N2, M] covariate4;
  matrix[N2, M] covariate5;
  matrix[N2, M] covariate6;
  int EpidemicStart[M]; // Currently always 31
  real pop[M]; // population number of countries
  real SI[N2]; // fixed pre-calculated serial intervals
}

transformed data {
  matrix[N2, M] matrix_0 = rep_matrix(0, N2, M);
}

parameters {
  // mortality
  real<lower=0> phi; // shape parameter for death distribution
  vector[M] ifr_std; // noise multiplying mean mortality
  real<lower=0,upper=3> log_ifr_mean;
  // R0 before interventions
  vector[M] mu_std;
  real<lower=2,upper=20> mu_mean;
  // covariate coefficients
  vector<lower = 0>[6] alpha_hier; // sudo parameter for the hier term for alpha
  // seed infections
  vector<lower = 0>[M] y_unit; // distribution of initial number of infecteds
  real<lower = 0> tau_unit; // hyperparamter for initial number of infecteds
}

transformed parameters {
  vector[6] alpha = alpha_hier - (log(1.05) / 6.0); // coefficients of covariates
  matrix[N2, M] prediction = matrix_0; // predicted number of daily infections
  matrix[N2, M] E_deaths = matrix_0; // predicted number of daily deaths
  // R0 with interventions (will be rescaled by initial R0 later)
  matrix[N2, M] Rt = exp(-(covariate1 * alpha[1]
                           + covariate2 * alpha[2]
                           + covariate3 * alpha[3]
                           + covariate4 * alpha[4]
                           + covariate5 * alpha[5]
                           + covariate6 * alpha[6]));
  real ifr_mean = 10^log_ifr_mean;
  real<lower = 0> tau = tau_unit / 0.03 * ifr_mean;
  vector<lower = 0>[M] y = tau * y_unit; // y ~ exponential(1/tau)
  vector<lower = 0>[M] ifr = exp(ifr_std * 0.1) * ifr_mean;
  vector<lower = 0>[M] mu = exp(mu_std * 0.1) * mu_mean;
  matrix[N2, M] cumm_sum = matrix_0; // total number of immune individuals
  
  for (m in 1:M) {
    for (i in 2:N0) {
      cumm_sum[i, m] = cumm_sum[i-1, m] + y[m]; 
    }
    // The loop below is faster than the code
    // prediction[1:N0, m] = rep_vector(y[m], N0);
    // because it doesn't require a new vector copy 
    // and loops without autodiff are very fast in Stan (Bob Carpenter)
    for (i in 1:N0) {
      prediction[i, m] = y[m]; // learn the number of cases in the first N0 days
    }
    Rt[, m] *= mu[m];
    
    for (i in (N0+1):N2) {
      real convolution=0;
      for(j in 1:(i-1)) {
        convolution += Rt[j,m] * prediction[j, m] * SI[i-j] *
                         ((pop[m]-cumm_sum[i-1,m]) / pop[m]);
      }
      cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i-1,m];
      prediction[i, m] = convolution;
    }
    
    E_deaths[1, m]= 1e-15 * prediction[1, m]; // zero really
    for (i in 2:N2){
      for(j in 1:(i-1)) {
        E_deaths[i, m] += prediction[j ,m] * f[i-j, m] / ifr[m];
      }
    }
  }
}

model {
  tau_unit ~ exponential(1);
  y_unit ~ exponential(1);
  ifr_std ~ std_normal();
  mu_std ~ std_normal();
  alpha_hier ~ gamma(.1667, 1);
  phi ~ normal(0, 5);
  for(m in 1:M){
    deaths[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:N[m], m], phi);
  }
}

generated quantities {
  // predictions in the absence of interventions
  // code similar to above just with interventions removed and 
  // 0 added to variable names
  
  matrix[N2, M] prediction0 = matrix_0; // predicted daily infections
  matrix[N2, M] E_deaths0  = matrix_0; // predicted daily deaths
  matrix[N2, M] Rt_adj = Rt; // Rt adjusted by immunity
  
  {
    matrix[N2,M] cumm_sum0 = matrix_0;
    for (m in 1:M){
      for (i in 2:N0){
        cumm_sum0[i, m] = cumm_sum0[i-1, m] + y[m]; 
      }
      prediction0[1:N0,m] = rep_vector(y[m], N0); 
      for (i in (N0+1):N2) {
        real convolution0 = 0;
        for(j in 1:(i-1)) {
          convolution0 += prediction0[j, m] * SI[i-j]; 
        }
        cumm_sum0[i,m] = cumm_sum0[i-1,m] + prediction0[i-1,m];
        prediction0[i, m] =  ((pop[m]-cumm_sum0[i,m]) / pop[m]) * mu[m] * convolution0;
        
        Rt_adj[i,m] = ((pop[m]-cumm_sum[i,m]) / pop[m]) * Rt[i,m];
      }
      
      E_deaths0[1, m] = uniform_rng(1e-16, 1e-15);
      for (i in 2:N2){
        for(j in 1:(i-1)){
          E_deaths0[i,m] += prediction0[j,m] * f[i-j,m] / ifr[m];
        }
      }
    }
  }
}
