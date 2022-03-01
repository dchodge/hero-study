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
    
    vector[n] ab_measure;
    vector[n] Time; 
}
transformed data{
    vector[n] ab_measure_t;
    ab_measure_t = ab_measure/10.0; 
}
parameters {
  vector[num_age] z_a;
  vector[num_eth] z_e;
  vector[num_gen] z_g;
  vector[num_gen] z_s;

  real<lower = 0>sigma_a_1;
  real<lower = 0>sigma_g_1;
  real<lower = 0>sigma_e_1;
  real<lower = 0>sigma_s_1;
  
  real<lower = 0> a;
  real<lower = 0> sigma;
}
transformed parameters{
  // regression model
  vector[n] mu;
  mu = a + z_s[symp]*sigma_s_1 + z_a[age_group]*sigma_a_1 + z_e[eth]*sigma_e_1 + z_g[gen]*sigma_g_1;
}
model {
  // likelihood
  ab_measure_t ~ normal(mu, sigma);
  // priors
  
  a ~ normal(0.6, 0.05);
  
  sigma ~ exponential(1);

  z_a[age_group] ~ normal(0, 1);
  sigma_a_1 ~ exponential(1);
  z_e[eth] ~ normal(0, 1);
  sigma_e_1 ~ exponential(1);
  z_g[gen] ~ normal(0, 1);
  sigma_g_1 ~ exponential(1);
  z_s[symp] ~ normal(0, 1);
  sigma_s_1 ~ exponential(1);
}
generated quantities{
  real marginal_i;
  vector[num_eth] marginal_e_i;
  vector[num_gen] marginal_g_i;
  vector[num_age] marginal_a_i;
  vector[num_symp] marginal_s_i;
  
  marginal_i = 0;
  marginal_e_i = rep_vector(0, num_eth);
  marginal_a_i = rep_vector(0, num_age);
  marginal_g_i = rep_vector(0, num_gen);
  marginal_s_i = rep_vector(0, num_symp);

  for (ai in 1:num_age){
    for (g in 1:num_gen){
      for (e in 1:num_eth){
        for (s in 1:num_symp){
          marginal_i += (a + z_s[s]*sigma_s_1 + z_a[ai]*sigma_a_1 + z_e[e]*sigma_e_1 + z_g[g]*sigma_g_1)*p_age[ai]*p_gen[g]*p_symp[s]*p_eth[e]*10;
          marginal_e_i[e] += (a + z_s[s]*sigma_s_1 + z_a[ai]*sigma_a_1 + z_e[e]*sigma_e_1 + z_g[g]*sigma_g_1)*p_age[ai]*p_gen[g]*p_symp[s]*10;
          marginal_g_i[g] += (a + z_s[s]*sigma_s_1 + z_a[ai]*sigma_a_1 + z_e[e]*sigma_e_1 + z_g[g]*sigma_g_1)*p_age[ai]*p_eth[e]*p_symp[s]*10;
          marginal_a_i[ai] += (a + z_s[s]*sigma_s_1 + z_a[ai]*sigma_a_1 + z_e[e]*sigma_e_1 + z_g[g]*sigma_g_1)*p_eth[e]*p_gen[g]*p_symp[s]*10;
          marginal_s_i[s] += (a + z_s[s]*sigma_s_1 + z_a[ai]*sigma_a_1 + z_e[e]*sigma_e_1 + z_g[g]*sigma_g_1)*p_eth[e]*p_gen[g]*p_age[ai]*10;
        }
      }
    }
  }
}

