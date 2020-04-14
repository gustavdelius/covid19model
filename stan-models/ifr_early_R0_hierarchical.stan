data {
  int <lower=1> M; // number of countries
  int <lower=1> N0; // number of days for which to impute infections
  int<lower=1> N[M]; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  int cases[N2,M]; // reported cases
  int deaths[N2, M]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  matrix[N2, M] f; // h * s * ifr
  matrix[N2, M] covariate1;
  matrix[N2, M] covariate2;
  matrix[N2, M] covariate3;
  matrix[N2, M] covariate4;
  matrix[N2, M] covariate5;
  matrix[N2, M] covariate6;
  int EpidemicStart[M];
  real pop[M];
  real SI[N2]; // fixed pre-calculated SI using emprical data from Neil
}

parameters {
  real<lower=0> alpha_hier[6]; // sudo parameter for the hier term for alpha
  real<lower=0> y_raw[M];
  real<lower=0> phi;
  real<lower=0> tau_unit;
  real<lower=0,upper=3> log_ifr_mean;
  real ifr_std[M];
  real<lower=2,upper=20> mu_mean;
  real mu_std[M];
}

transformed parameters {
    real alpha[6];
    matrix[N2, M] prediction = rep_matrix(0,N2,M);
    matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
    matrix[N2, M] Rt = rep_matrix(0,N2,M);
    real y[M];
    real ifr[M];
    real mu[M];
    
    real ifr_mean = 10^log_ifr_mean;
    real<lower = 0> tau = tau_unit / 0.03 * ifr_mean;

      matrix[N2,M] cumm_sum = rep_matrix(0,N2,M);
      for(i in 1:6){
        alpha[i] = alpha_hier[i] - ( log(1.05) / 6.0 );
      }
      for (m in 1:M){
        y[m] = tau * y_raw[m];
        ifr[m] = exp(ifr_std[m] * 0.1) * ifr_mean;
        mu[m] = exp(mu_std[m] * 0.1) * mu_mean;
        for (i in 2:N0){
          cumm_sum[i,m] = cumm_sum[i-1,m] + y[m]; 
        }
        prediction[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
        
        Rt[,m] = mu[m] * exp( covariate1[,m] * (-alpha[1]) + covariate2[,m] * (-alpha[2]) +
          covariate3[,m] * (-alpha[3]) + covariate4[,m] * (-alpha[4]) + covariate5[,m] * (-alpha[5]) + 
          covariate6[,m] * (-alpha[6]) );
        for (i in (N0+1):N2) {
          real convolution=0;
          for(j in 1:(i-1)) {
            convolution += Rt[j,m] * prediction[j, m] * SI[i-j] *
              ((pop[m]-cumm_sum[i-1,m]) / pop[m]);
          }
          cumm_sum[i,m] = cumm_sum[i-1,m] + prediction[i-1,m];
          prediction[i, m] = convolution;
        }
        
        E_deaths[1, m]= 1e-15 * prediction[1,m];
        for (i in 2:N2){
          for(j in 1:(i-1)){
            E_deaths[i,m] += prediction[j,m] * f[i-j,m] / ifr[m] ;
          }
        }
      }
}
model {
  tau_unit ~ exponential(1);
  y_raw ~ exponential(1);
  ifr_std ~ std_normal();
  mu_std ~ std_normal();
  phi ~ normal(0,5);
  alpha_hier ~ gamma(.1667,1);
  for(m in 1:M){
    deaths[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:N[m], m], phi);
   }
}

generated quantities {
    matrix[N2, M] prediction0 = rep_matrix(0,N2,M);
    matrix[N2, M] E_deaths0  = rep_matrix(0,N2,M);
    matrix[N2, M] Rt_adj = Rt;
    
    {
      matrix[N2,M] cumm_sum0 = rep_matrix(0,N2,M);
      for (m in 1:M){
         for (i in 2:N0){
          cumm_sum0[i,m] = cumm_sum0[i-1,m] + y[m]; 
        }
        prediction0[1:N0,m] = rep_vector(y[m],N0); 
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
