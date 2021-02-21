data
{
  int n;
  int n_cov;
  int y_se;
  int n_se;
  int y_sp;
  int n_sp;
  int response[n];
  int tot[n];
  int covariates[n, n_cov];
  int cov_size[n_cov + 1];
  int cov_size_tot;
  int indicies[prod(cov_size[2:]), n_cov];
  matrix[max(cov_size), n_cov] prop;
}
transformed data{
  int l;

  l = prod(cov_size[2:]);
}
parameters
{
  real beta_0;
  vector[cov_size_tot] beta;
  vector<lower = 0>[n_cov] sigma;
  real se;
  real sp;
}
transformed parameters
{
  real<lower = 0, upper = 1> se_t; 
  real<lower = 0, upper = 1> sp_t;
  se_t = inv_logit(se);
  sp_t = inv_logit(sp);
}
model
{
  vector[n] pi_true;
  vector[n] pi_obs;
  int ii[n_cov];

  for (i in 1:n){
    for (j in 1:n_cov){
      ii[j] = cov_size[j] + covariates[i][j];
    }
    pi_true[i] = inv_logit(beta_0 + sum(beta[ii].*sigma));
  }
  
  beta_0 ~ normal(0, 1);
  pi_obs = pi_true*se_t + (1-pi_true)*(1-sp_t);
  response ~ binomial(tot, pi_obs);
  y_se ~ binomial(n_se, se_t);
  y_sp ~ binomial(n_sp, sp_t);
  beta ~ normal(0, 1);
  sigma ~ exponential(1);

  se ~ normal(5, 1);
  sp ~ normal(5, 1);
}
generated quantities
{
  real prod_val;
  matrix[max(cov_size), n_cov] marginal;
  
  int index[n_cov];
  

  for (c1 in 1:max(cov_size)){
     for (c2 in 1:n_cov){
      marginal[c1, c2] = 0;
     }
  }

  for (c in 1:n_cov){
    for (ii in 1:l){
      prod_val =  1;
      for (c2 in 1:n_cov){
        prod_val *= prop[indicies[ii, c2], c2];
      }
      for (j in 1:n_cov){
        index[j] = cov_size[j] + indicies[ii][j];
      }
      marginal[indicies[ii, c], c] += inv_logit(beta_0 + sum(beta[index].*sigma))*prod_val;
    }
  }
  
  for (c1 in 1:max(cov_size)){
     for (c2 in 1:n_cov){
      marginal[c1, c2] = marginal[c1, c2]/prop[c1, c2];
     }
  }
}
