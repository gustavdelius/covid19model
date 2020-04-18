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
  matrix[N2, M] matrix_1 = rep_matrix(1, N2, M);
}

parameters {
  // mortality
  real<lower=0> phi; // shape parameter for death distribution
  vector<lower=0>[M] ifr_noise; // noise multiplying mean mortality
  real log_infecteds_multiplier; // We will _divide_ ifr by exp(log_infecteds_multiplier)
  // R0 before interventions
  vector<lower=0>[M] mu; // initial value of R0
  real<lower=0> kappa; // standard deviation of initial R0
  real log_mu_factor;
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
                           + covariate6 * alpha[6])
                         + log_mu_factor);
  // The IFR calculated by Imperial group will be multiplied by the random ifr_factor
  vector<lower = 0>[M] ifr_factor = ifr_noise / exp(log_infecteds_multiplier); 
  // Number of initial infecteds should also be multiplied by infectes_multiplier,
  // hence tau ~ exponential(0.03) * infecteds_multiplier
  real<lower = 0> tau = tau_unit / 0.03 * exp(log_infecteds_multiplier); 
  vector<lower = 0>[M] y = tau * y_unit; // y ~ exponential(1/tau)
  matrix[N2, M] susceptible = matrix_1; // proportion of susceptibles in population
  
  for (m in 1:M) {
    prediction[1, m] = y[m];
    for (i in 2:N0) {
      prediction[i, m] = y[m]; // keeping same number of cases in the first N0 days
      susceptible[i, m] = susceptible[i-1, m] - y[m] / pop[m];
    }
    Rt[, m] *= mu[m];
    
    for (i in (N0 + 1):N2) {
      real convolution = 0;
      for(j in 1:(i-1)) {
        convolution += prediction[j, m] * SI[i-j];
      }
      susceptible[i, m] = susceptible[i-1, m] - prediction[i-1, m] / pop[m];
      prediction[i, m] = susceptible[i, m] * Rt[i, m] * convolution;
    }
    E_deaths[1, m]= 1e-15 * prediction[1, m]; // zero really
    for (i in 2:N2) {
      for(j in 1:(i-1)) {
        E_deaths[i, m] += prediction[j, m] * f[i-j, m] * ifr_factor[m];
      }
    }
  }
}

model {
  tau_unit ~ exponential(1); // tau ~ exponential(0.03)
  y_unit ~ exponential(1); // y ~ exponential(1/tau)
  kappa ~ normal(0, 0.5);
  mu ~ normal(3.28, kappa); // citation: https://academic.oup.com/jtm/article/27/2/taaa021/5735319
  log_mu_factor ~ std_normal();
  alpha_hier ~ gamma(.1667, 1);
  ifr_noise ~ normal(1, 0.1);
  log_infecteds_multiplier ~ normal(0, 2);
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
  matrix[N2, M] Rt_adj = Rt .* susceptible; // Effective R including immunity
  
  {
    matrix[N2, M] susceptible0 = matrix_1; // proportion of susceptibles in population
    
    for (m in 1:M) {
      prediction0[1, m] = y[m];
      for (i in 2:N0) {
        prediction0[i, m] = y[m]; // keeping same number of cases in the first N0 days
        susceptible0[i, m] = susceptible0[i-1, m] - y[m] / pop[m];
      }
      
      for (i in (N0 + 1):N2) {
        real convolution = 0;
        for(j in 1:(i-1)) {
          convolution += mu[m] * prediction0[j, m] * SI[i-j];
        }
        susceptible0[i, m] = susceptible0[i-1, m] - prediction0[i-1, m] / pop[m];
        prediction0[i, m] = susceptible0[i, m] * convolution;
      }
      
      E_deaths0[1, m] = uniform_rng(1e-16, 1e-15);
      for (i in 2:N2){
        for(j in 1:(i-1)){
          E_deaths0[i, m] += prediction0[j, m] * f[i-j, m] * ifr_factor[m];
        }
      }
    }
  }
}

