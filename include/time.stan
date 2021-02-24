data { 
  int num_data;
  vector[100] time_seq;
  vector[num_data] time; 
  vector[num_data] measure;
} 
transformed data {
  vector[num_data] time_t; 
  vector[num_data] measure_t;
  real threshold_measure;

  time_t = time/100.0;
  measure_t = measure/10.0;
  threshold_measure = min(measure);
}
parameters { 
  // linear model
  real a;
  real b;
  real<lower= 0> sigma;
} 
model { 
  measure_t ~ normal(a + b*time_t, sigma); 
  a ~ normal(0.5, 0.2);
  b ~ normal(0, 1);
  sigma ~ exponential(1);
}
generated quantities{
  vector[num_data] log_lik;

  vector[100] mean_pred;
  real pred[100];
  
  real time_revert_raw;
  real<lower = -1, upper = 365*10 + 1> time_revert;

  real t0;
  real y0;
  real t_half;
  
  t0 = 14.0/100.0;
  y0 = a + b*t0;
  
  mean_pred = (a + b*time_seq/100)*10;
  pred = normal_rng(mean_pred, sigma);
  
  time_revert_raw = (threshold_measure/10.0 - a)*100/b;
  if (time_revert_raw < 0)
    time_revert = -1;
  else if (time_revert_raw > 365*10)
    time_revert = 365*10 + 1;
  else 
    time_revert = time_revert_raw;
    
  for (n in 1:num_data) {
    log_lik[n] = normal_lpdf(measure_t[n] | a + b*time_t[n], sigma);
  }
  
  t_half = -log(2)*10/b;
}
