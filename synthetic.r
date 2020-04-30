library(rstan)
load("results/base-stanfit.Rdata")

# data
M <- stan_data$M  # number of countries
N <- stan_data$N  # days of observed data for country m. each entry must be <= N2
EpidemicStart <- stan_data$EpidemicStart

# parameters
# We choose those from the posterior sample with the highest density
out <- rstan::extract(fit)
idx <- which.max(out$lp__)
phi <- out$phi[idx]
E_deaths <- out$E_deaths[idx, , ]

# model
for(m in 1:M){
  for (i in EpidemicStart[m]:N[m]) {
    stan_data$deaths[i, m] <- rnbinom(1, mu = E_deaths[i, m], size = phi)
    }
}

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)
StanModel <- "base_ifr"
m = stan_model(paste0('stan-models/',StanModel,'.stan'))
fit = sampling(
  m,
  data = stan_data,
  iter = 1800,
  warmup = 1000,
  chains = 5,
  thin = 1,
  control = list(adapt_delta = 0.95, max_treedepth = 15)
)