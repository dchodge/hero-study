data {
    int n;
    int num_age;
    int num_eth;
    int num_gen;
    int num_symp;

    real p_gen[num_gen];
    real p_age[num_age];
    real p_eth[num_eth];
    real p_symp[num_symp];

    int age_group[n]; 
    int eth[n]; 
    int gen[n]; 
    int symp[n];
    
    real ab_measure[n]; 
    vector[n] Time; 
}
transformed data{
  vector[n] Time_t; 

  Time_t = Time/28.0;
}
parameters {
  vector[num_age] x_a;
  vector[num_eth] x_e;
  vector[num_gen] x_g;
  vector[num_symp] x_s;

  real<lower = 0>sigma;
  real<lower = 0>sigma_a_2;
  real<lower = 0>sigma_g_2;
  real<lower = 0>sigma_e_2;
  real<lower = 0>sigma_s_2;

  real b;
}
transformed parameters{
  // regression model
  vector[n] mu;
  mu = (Time_t).*(b + x_s[symp]*sigma_s_2 + x_a[age_group]*sigma_a_2 + x_e[eth]*sigma_e_2 + x_g[gen]*sigma_g_2);

}
model {
  // likelihood
  ab_measure ~ normal(mu, sigma);
  // priors
  b + mean(x_s)*sigma_s_2 + mean(x_a)*sigma_a_2 + mean(x_g)*sigma_g_2 + mean(x_e)*sigma_e_2 ~ normal(0, 0.5);

  sigma ~ exponential(1);
  
  x_a[age_group] ~ normal(0, 1);
  sigma_a_2 ~ exponential(1);
  x_e[eth] ~ normal(0, 1);
  sigma_e_2 ~ exponential(1);
  x_g[gen] ~ normal(0, 1);
  sigma_g_2 ~ exponential(1);
  x_s[symp] ~ normal(0, 1);
  sigma_s_2 ~ exponential(1);
}
generated quantities{
// get out the covariance and correlation matrices

  real marginal_s;
  vector[num_eth] marginal_e_s;
  vector[num_gen] marginal_g_s;
  vector[num_age] marginal_a_s;
  vector[num_symp] marginal_s_s;


  marginal_s = 0;
  marginal_e_s = rep_vector(0, num_eth);
  marginal_a_s = rep_vector(0, num_age);
  marginal_g_s = rep_vector(0, num_gen);
  marginal_s_s = rep_vector(0, num_symp);

  for (ai in 1:num_age){
    for (g in 1:num_gen){
      for (e in 1:num_eth){
        for (s in 1:num_symp){
          marginal_s += (b + x_s[s]*sigma_s_2 + x_a[ai]*sigma_a_2 + x_e[e]*sigma_e_2 + x_g[g]*sigma_g_2)*p_age[ai]*p_gen[g]*p_symp[s]*p_eth[e];
          marginal_e_s[e] += (b + x_s[s]*sigma_s_2 + x_a[ai]*sigma_a_2 + x_e[e]*sigma_e_2 + x_g[g]*sigma_g_2)*p_age[ai]*p_gen[g]*p_symp[s];
          marginal_g_s[g] += (b + x_s[s]*sigma_s_2 + x_a[ai]*sigma_a_2 + x_e[e]*sigma_e_2 + x_g[g]*sigma_g_2)*p_age[ai]*p_eth[e]*p_symp[s];
          marginal_a_s[ai] += (b + x_s[s]*sigma_s_2 + x_a[ai]*sigma_a_2 + x_e[e]*sigma_e_2 + x_g[g]*sigma_g_2)*p_eth[e]*p_gen[g]*p_symp[s];
          marginal_s_s[s] += (b + x_s[s]*sigma_s_2 + x_a[ai]*sigma_a_2 + x_e[e]*sigma_e_2 + x_g[g]*sigma_g_2)*p_eth[e]*p_gen[g]*p_age[ai];
        }
      }
    }
  }
}
