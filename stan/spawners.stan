data {
  int                    n_priors;        // n priors
  matrix[n_priors,2]     priors;          // rows: priors; cols: p1, p2
  int                    n_years;         // n years
  int                    n_obs;          // number of spawner count obs
  
  int<lower=1, upper=n_years>     year[n_obs];
  int                    day[n_obs];
  int                    live_counts[n_obs]; //live spawner counts
}
parameters {
  // group distributions (years)
  real timing_mu;         
  real<lower=0> timing_sigma;           // variation across years
  
  real      spread_arrive_mu;     // stdev entry
  real<lower=0> spread_arrive_sigma;     
  
  real      spread_death_mu;      // stdev death
  real<lower=0> spread_death_sigma;     
  
  // simple
  real     residence_raw;  // lag entry to exit
  real<lower=0>         phi_live;   // dispersion parameter on live fish sampling
  
  // instance in group, multi-level
  vector[n_years] log_run;
  vector[n_years] timing_z;
  vector[n_years] spread_arrive_z;
  vector[n_years] spread_death_z;
}

transformed parameters{
  //non-centered priors
  vector[n_years] timing;
  vector[n_years] spread_arrive;
  vector[n_years] spread_death;
  
  
  for (y in 1:n_years){
    timing[y] = timing_mu + timing_sigma * timing_z[y];
    spread_arrive[y] = exp(spread_arrive_mu + spread_arrive_sigma * spread_arrive_z[y]) + 1;
    spread_death[y] = exp(spread_death_mu + spread_death_sigma * spread_death_z[y]) + 1;
  }
  
  real residence;
  residence = exp(residence_raw) + 1;
  
  //process model
  // live
  vector[n_obs] mu_live; // predictions from model
  for(i in 1:n_obs){
    real mean_arrival = timing[year[i]];
    real entered = normal_cdf(day[i], mean_arrival, spread_arrive[year[i]]);  
    real exited  = normal_cdf(day[i], mean_arrival + residence, spread_death[year[i]]);
    mu_live[i] = exp(log_run[year[i]] +log1p(entered * (1 - exited) * 1e4) - log(1e4));
  }
}

model {
  // local variables
  // group prior PDDs
  
  timing_mu ~ normal(priors[2,1], priors[2,2]); 
  timing_sigma ~ exponential(priors[3,2]);                
  timing_z ~ normal(0, 1);             
  
  spread_arrive_mu  ~ normal(priors[4,1], priors[4,2]);  // 
  spread_arrive_sigma ~ exponential(priors[5,2]);
  spread_arrive_z ~ normal(0, 1);
  
  spread_death_mu   ~ normal(priors[4,1], priors[4,2]);  // 
  spread_death_sigma ~ exponential(priors[5,2]);
  spread_death_z ~ normal(0, 1);
  
  // simple PDD
  log_run ~ normal(priors[1,1], priors[1,2]);  // mean  run log
  residence_raw  ~ normal(priors[6,1], priors[6,2]);  // same all years
  log(phi_live)   ~ normal(priors[7,1], priors[7,2]);
  
  
  //likelihood
  target += neg_binomial_2_lpmf(live_counts | mu_live, phi_live);
}

generated quantities{
  vector[n_obs] log_lik;
  for(i in 1:n_obs)
  log_lik[i] = neg_binomial_2_lpmf(live_counts[i] | mu_live[i], phi_live);
}
