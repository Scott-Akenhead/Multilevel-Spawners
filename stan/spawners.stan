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
  real arrival_mu;         
  real<lower=0> arrival_sigma;           // variation across years
  
  real      arrival_spread_mu;     // stdev entry
  real<lower=0> arrival_spread_sigma;     
  
  real      death_spread_mu;      // stdev death
  real<lower=0> death_spread_sigma;     
  
  // simple
  real     residence_raw;  // lag entry to exit
  real<lower=0>         live_phi;   // dispersion parameter on live fish sampling
  
  // instance in group, multi-level
  vector[n_years] log_run;
  vector[n_years] arrival_z;
  vector[n_years] arrival_spread_z;
  vector[n_years] death_spread_z;
}

transformed parameters{
  //non-centered priors
  vector[n_years] arrival;
  vector[n_years] arrival_spread;
  vector[n_years] death_spread;
  
  
  for (y in 1:n_years){
    arrival[y] = arrival_mu + arrival_sigma * arrival_z[y];
    arrival_spread[y] = exp(arrival_spread_mu + arrival_spread_sigma * arrival_spread_z[y]) + 1;
    death_spread[y] = exp(death_spread_mu + death_spread_sigma * death_spread_z[y]) + 1;
  }
  
  real residence;
  residence = exp(residence_raw) + 1;
  
  //process model
  // live
  vector[n_obs] live_mu; // predictions from model
  for(i in 1:n_obs){
    real mean_arrival = arrival[year[i]];
    real entered = normal_cdf(day[i], mean_arrival, arrival_spread[year[i]]);  
    real exited  = normal_cdf(day[i], mean_arrival + residence, death_spread[year[i]]);
    live_mu[i] = exp(log_run[year[i]] +log1p(entered * (1 - exited) * 1e4) - log(1e4));
  }
}

model {
  // local variables
  // group prior PDDs
  
  arrival_mu ~ normal(priors[2,1], priors[2,2]); 
  arrival_sigma ~ exponential(priors[3,2]);                
  arrival_z ~ normal(0, 1);             
  
  arrival_spread_mu  ~ normal(priors[4,1], priors[4,2]);  // 
  arrival_spread_sigma ~ exponential(priors[5,2]);
  arrival_spread_z ~ normal(0, 1);
  
  death_spread_mu   ~ normal(priors[4,1], priors[4,2]);  // 
  death_spread_sigma ~ exponential(priors[5,2]);
  death_spread_z ~ normal(0, 1);
  
  // simple PDD
  log_run ~ normal(priors[1,1], priors[1,2]);  // mean  run log
  residence_raw  ~ normal(priors[6,1], priors[6,2]);  // same all years
  log(live_phi)   ~ normal(priors[7,1], priors[7,2]);
  
  
  //likelihood
  target += neg_binomial_2_lpmf(live_counts | live_mu, live_phi);
}

generated quantities{
  vector[n_obs] log_lik;
  for(i in 1:n_obs)
  log_lik[i] = neg_binomial_2_lpmf(live_counts[i] | live_mu[i], live_phi);
}
