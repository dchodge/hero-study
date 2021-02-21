# Run this to dave all the data for fitting.
save_all_datafit <- function() {
  data_full <- output_rdata()
  datafit_prev <- get_datafit_prev()
  datafit_asymp <- get_datafit_asymp()
  datafit_start <- get_datafit_start()
  datafit_change <- get_datafit_change()
  datafit_seropos <- get_datafit_seropos()
  datafit_seroneg <- get_datafit_seroneg()
  datafit_od_val <- get_datafit_od_val(0.175, 0.1905)
  datafit_control <- get_datafit_control()
  save(data_full, datafit_prev, datafit_asymp, datafit_start, datafit_change,
       datafit_seropos, datafit_seroneg, datafit_od_val, datafit_control, file = here::here("data", "datafit.RData"))
}

# fit all the models in the main regression (takes a couple of hours)
run_fit_mcmc_all <- function() {
  prev <- init_model_prev("prev", "sero_pos", datafit_prev, "sero_prev")
  asymp <- init_model_prev("asymp", "asymp_pos", datafit_asymp, "asymp")
  inter_start <- c("marginal_a_i", "marginal_e_i", "marginal_g_i", "marginal_s_i")
  slope_change <- c("marginal_a_s", "marginal_e_s", "marginal_g_s", "marginal_s_s")
  start <- init_model_ab_kin(name = "start",
                             stan_file = "start",
                             datafit = datafit_start,
                             marginal_vars = inter_start)
  change <- init_model_ab_kin(name = "change",
                              stan_file = "change",
                              datafit = datafit_change,
                              marginal_vars = slope_change)
  # Run the mcmc sampler on these fits
  fit(prev)
  fit(asymp)
  fit(start)
  fit(change)
}

# Load all the data for fitting
load_fit_mcmc_all <- function() {
  data(fit_mcmc_prev)
  data(fit_mcmc_asymp)
  data(fit_mcmc_start)
  data(fit_mcmc_change)
}

# Save all the data for plotting
save_all_dataplt <- function() {
  dataplt_prev <- get_dataplt_prev()
  dataplt_asymp <- get_dataplt_asymp()
  dataplt_start <- get_dataplt_start()
  dataplt_change <- get_dataplt_change()
  dataplt_revert <- get_dataplt_revert()
  dataplt_sens <- get_dataplt_sens(seq(70, 100, 1))
  dataplt_spec <- get_dataplt_spec(seq(70, 100, 1))
  save(dataplt_prev, dataplt_asymp, dataplt_start, dataplt_change, dataplt_revert,
       dataplt_sens, dataplt_spec, file = here::here("data", "dataplt.RData"))
}

get_mean_ci <- function(vec) {
  l1 <- vec %>% mean %>% setNames("Mean")
  l2 <- vec %>% quantile(c(0.025, 0.975))
  c(l1, l2)
}

mutate_quantile <- function(x, probs = c(0.025, 0.5, 0.975)) {
  tibble(value = quantile(x, probs), probs = probs)
}