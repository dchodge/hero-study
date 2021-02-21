## Plot Supplementary Figure 1
plot_dates <- function(height, width) {
  df_dates <- data.frame(
    "Symptom Onset" = as.Date(data_full$symp$inf_date, origin = "2020-01-01")[seq(1, 326, 2)],
    "First bleed" = as.Date(data_full$symp$date_of_sample, origin = "2020-01-01")[seq(1, 326, 2)],
    "Second bleed" = as.Date(data_full$symp$date_of_sample, origin = "2020-01-01")[seq(2, 326, 2)]) %>%
    tidyr::gather(event, value)
  ggplot(data = df_dates) +
    geom_histogram(aes(x = value, fill = event), color = "black") +
    labs(x = "Date (2020)", y = "Count", fill = "Event") + theme_bw()
  ggsave(filename = here::here("outputs", "figs", "supp_timings.pdf"), height = height, width = width)
}


############## Find correlations between ab measures ###############
get_df_correlation_model <- function() {
  df_ab_long_d <- datafit_start %>% select(IgG_S, IgG_N, IgA_S_Serum, IgA_N_Serum,
                                        age_group,  time)
  df_ab_long_d$sample_no <- purrr::map(1:(nrow(df_ab_long_d) / 2), ~rep(.x, 2)) %>% unlist
  df_ab_long_d$age_group <- as.numeric(df_ab_long_d$age_group)
  df_ab_long_d <- df_ab_long_d[seq(1, nrow(df_ab_long_d), 2), ]
  df_ab_long_d$IgG_S <- scale(df_ab_long_d$IgG_S)
  df_ab_long_d$IgG_N <- scale(df_ab_long_d$IgG_N)
  df_ab_long_d$IgA_S_Serum <- scale(df_ab_long_d$IgA_S_Serum)
  df_ab_long_d$IgA_N_Serum <- scale(df_ab_long_d$IgA_N_Serum)
  data_measure <- df_ab_long_d %>% select(IgG_S, IgG_N, IgA_S_Serum, IgA_N_Serum)
  data_measure
}

get_rsquared <- function() {
  data_measure <- get_df_correlation_model()
  measure_names <- c("IgG-S", "IgG-N", "IgA-serum-S", "IgA-serum-N")

  corr_val <- matrix(, 4, 4)
  data_r <- data.frame()
  for (i in 1:4) {
    for (j in 1:4) {
      data_r_t <- data_measure[, c(i, j)] %>% setNames(c("x", "y"))
      corr_val[i, j] <- round(cor(data_r_t[["x"]], data_r_t[["y"]])^2, 2)
      data_r_t$measure1 <- measure_names[i]
      data_r_t$measure2 <-  measure_names[j]
      data_r <- rbind(data_r, data_r_t)
    }
  }
  corr_val_df <- corr_val %>%
    as.data.frame %>%
    setNames(measure_names) %>%
    tidyr::gather(key = "measure1", value) %>%
    mutate(measure2 = rep(measure_names, 4))
  compare_ab <- list(corr = corr_val_df, data = data_r)
  compare_ab
}

############## Find correlations between ab measures ###############
plot_correlations <- function(height, width) {
  compare_ab <- get_rsquared()
  measure_names <- c("IgG-S", "IgG-N", "IgA-serum-S", "IgA-serum-N")
  ggplot(data = compare_ab$data, aes(x = x, y = y)) +
    geom_point(alpha = 0.6, size = 0.1) +
    facet_grid(rows = vars(factor(measure1, levels = measure_names)),
      cols = vars(factor(measure2, levels = measure_names))) +
    geom_label(data = compare_ab$corr, aes(x = 0, y = 2, label = value)) +
    labs(x = "Antibody measure", y = "Antibody measure")
  ggsave(filename = here::here("outputs", "figs", "supp_corr.pdf"), height = height, width = width)
}





stan_model_linear_fit <- "
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
"

fit_mcmc_time_pso <- function() {
  m4_time <- list()
  # fit IgG and IgA serum
  measure_list <- c("IgG_S_unit", "IgG_N_unit", "IgA_Serum_S_unit", "IgA_Serum_N_unit")
  time_seq <- c(14:113) # generating input
  for (measure_unit in measure_list) {
    data_measure <- datafit_abkin %>% select(days_since_inf, measure_unit) %>% na.omit
    time <- data_measure$days_since_inf
    num_data <- length(time)
    measure <- data_measure[[measure_unit]] # generating inputs
    data_list <- list(num_data = num_data, time_seq = time_seq, time = time, measure = measure)
    m4_time[[measure_unit]] <- stan(model_code = stan_model_linear_fit, data = data_list, chains = 4, cores = 4)
  }
  save(m4_time,  file = here::here("data", "fit_mcmc_time_pso.RData"))
}