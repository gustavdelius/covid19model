data {
  int <lower=1> M; // number of countries
  int<lower=1> Pop[M]; // population of countries
  int<lower=1> N0; // number of days for which to impute infections
  int<lower=1> N[M]; // days of observed data for country m. each entry must be <= N2
  int<lower=1> N2; // days of observed data + # of days to forecast
  real<lower=0> x[N2]; // index of days (starting at 1)
  int cases[N2,M]; // reported cases
  int deaths[N2, M]; // reported deaths -- the rows with i > N contain -1 and should be ignored
  real F[N2]; // Distribution of time to death
  matrix[N2, M] covariate1;
  matrix[N2, M] covariate2;
  matrix[N2, M] covariate3;
  matrix[N2, M] covariate4;
  matrix[N2, M] covariate5;
  matrix[N2, M] covariate6;
  int EpidemicStart[M];
  real SI[N2]; // fixed pre-calculated SI using emprical data from Neil
}

parameters {
  real<lower=0> mu[M]; // intercept for Rt
  real<lower=0> alpha[6]; // the hier term
  real<lower=0> kappa;
  vector<lower=0>[M] y_raw;
  real<lower=0> phi;
  real<lower=0> tau;
  vector<lower=1, upper=9>[M] ifr_log; // log10 of infection-fatality ratio
}

transformed parameters {
    vector<lower = 0>[M] y = tau * y_raw ;
    vector<upper = 1>[M] fir = rep_vector(0, M);
    matrix[N2, M] f = rep_matrix(0,N2,M);
    matrix[N2, M] tc = rep_matrix(0,N2,M);
    matrix[N2, M] prediction = rep_matrix(0,N2,M);
    matrix[N2, M] E_deaths  = rep_matrix(0,N2,M);
    matrix[N2, M] Rt = rep_matrix(0,N2,M);
    for (m in 1:M) {
        vector[N2] h = rep_vector(1, N2);
        vector[N2] s = rep_vector(1, N2);
        fir[m] = 10^(-ifr_log[m]);
        y[m] = tau * y_raw[m] / fir[m];
        h[1] = fir[m] * F[1];
        s[1] = 1;
        f[1, m] = h[1];
        for(i in 2:N2) {
            h[i] = fir[m] * (F[i] - F[i-1]) / (1 - fir[m] * F[i]);
            s[i] = s[i-1] * (1 - h[i-1]);
            f[i, m] = s[i] * h[i];
        }
        prediction[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
        Rt[,m] = mu[m] * exp(covariate1[,m] * (-alpha[1]) + covariate2[,m] * (-alpha[2]) +
        covariate3[,m] * (-alpha[3])+ covariate4[,m] * (-alpha[4]) + covariate5[,m] * (-alpha[5]) + 
            covariate6[,m] * (-alpha[6])); // + GP[i]); // to_vector(x) * time_effect
        for (i in (N0+1):N2) {
            real convolution=0;
            for(j in 1:(i-1)) {
                convolution += prediction[j, m]*SI[i-j]; // Correctd 22nd March
            }
            tc[i, m] = tc[i-1, m] + prediction[i-1, m];
            prediction[i, m] = Rt[i,m] * convolution * (1 - tc[i, m] / Pop[m]);
        }
      
        E_deaths[1, m]= 1e-9;
        for (i in 2:N2){
            E_deaths[i,m]= 0;
            for(j in 1:(i-1)){
                E_deaths[i,m] += prediction[j,m]*f[i-j,m];
            }
        }
    }
   /* for(m in 1:M) {
     for(i in 1:N[m]) {
      LowerBound[i,m] = prediction[i,m] * 10 - cases[i,m];
     }
    }*/

}
model {
  tau ~ exponential(0.03);
  y_raw ~ exponential(1);  // implies y ~ exponential(1 / tau);
  phi ~ normal(0,5);
  kappa ~ normal(0,0.5);
  mu ~ normal(2.4, kappa); // citation needed 
  alpha ~ gamma(.5,1);
  for(m in 1:M){
    deaths[EpidemicStart[m]:N[m], m] ~ neg_binomial_2(E_deaths[EpidemicStart[m]:N[m], m], phi);
  }
}

generated quantities {
    matrix[N2, M] lp0 = rep_matrix(1000,N2,M); // log-probability for LOO for the counterfactual model
    matrix[N2, M] lp1 = rep_matrix(1000,N2,M); // log-probability for LOO for the main model
    matrix[N2, M] prediction0 = rep_matrix(0,N2,M);
    matrix[N2, M] E_deaths0  = rep_matrix(0,N2,M);
    for (m in 1:M){
      prediction0[1:N0,m] = rep_vector(y[m],N0); // learn the number of cases in the first N0 days
      for (i in (N0+1):N2) {
        real convolution0=0;
        for(j in 1:(i-1)) {
          convolution0 += prediction0[j, m]*SI[i-j]; // Correctd 22nd March
        }
        prediction0[i, m] = mu[m] * convolution0;
      }
      
      E_deaths0[1, m]= 1e-9;
      for (i in 2:N2){
        E_deaths0[i,m]= 0;
        for(j in 1:(i-1)){
          E_deaths0[i,m] += prediction0[j,m]*f[i-j,m];
        }
      }
      for(i in 1:N[m]){
        lp0[i,m] = neg_binomial_2_lpmf(deaths[i,m] | E_deaths[i,m],phi); 
        lp1[i,m] = neg_binomial_2_lpmf(deaths[i,m] | E_deaths0[i,m],phi); 
      }
    }

}
