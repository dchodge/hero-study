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

#### Fit the time model

fit_time_supp <- function() {

  m_time <- stan_model(file = here::here("include", "time.stan"))
  m4_time <- list()
  # fit IgG and IgA serum
  measure_list <- c("IgG_S", "IgG_N", "IgA_S_Serum", "IgA_N_Serum")
  time_seq <- c(14:113) # generating input
  for (measure_unit in measure_list) {
    data_measure <- get_datafit_time() %>% select(days_since_inf, measure_unit) %>% na.omit
    time_pso <- data_measure$days_since_inf
    num_data <- length(time_pso)
    measure <- data_measure[[measure_unit]] # generating inputs
    data_list <- list(num_data = num_data, time_seq = time_seq, time = time_pso, measure = measure)
    m4_time[[measure_unit]] <- sampling(m_time, data = data_list, chains = 4, cores = 4,
        sample_file = here::here(paste0("include/time/", measure_list[i])),
        diagnostic_file = here::here(paste0("include/time/", measure_list[i], "_diag"))
    )
  }
  save(m4_time,  file = here::here("data", "fit_mcmc_time.RData"))
}

plot_time_supp <- function(height, width) {
  if (!exists("m4_time")) {
    data(fit_mcmc_time)
  }

  measure_list <- c("IgG_S", "IgG_N", "IgA_S_Serum", "IgA_N_Serum")
  longdata <- (1:4 %>% map(~ (tidybayes::spread_draws(m4_time[[.x]], a, b) %>%
      mutate(antibody = measure_list[.x])))) %>%
    bind_rows
  longdatamean <- longdata %>% group_by(antibody) %>% summarise(a = mean(a), b = mean(b))

  data_time_plot <- get_datafit_time() %>%
    pivot_longer(cols = all_of(measure_list), names_to = "antibody", values_to =  "antibody_measure")

  longdata %>%
    ggplot() +
    geom_abline(aes(intercept = a * 10, slope = b / 10), alpha = 0.05, size = 0.05, color = "red") +
    geom_abline(data = longdatamean,
      aes(intercept = a * 10, slope = b / 10), color = "darkred", size = 2) +
    geom_point(data = data_time_plot, aes(x = days_since_inf, y = antibody_measure), shape = 5, alpha = 0.8) +
    facet_wrap(vars(factor(antibody, levels = measure_list)), nrow = 2) + theme_bw() + theme(aspect.ratio = 0.7) +
    labs(x = "Days post-symptom onset", y = "Antibody measure (logAU)")
  ggsave(filename = here::here("outputs", "figs", "supp_time_pso.pdf"), height = height, width = width)
}


plot_od_sens_curve <- function(height, width) {
  df_label <- data.frame(
    x = c(50, 50, 50, 50),
    y = c(0.3, 2.1, 0.3, 3),
    protein = c("S", "S", "N", "N"))
  datafit_od_curve <- get_datafit_control()
  ggplot(data = datafit_od_curve, aes(y = od)) +
    geom_point(aes(x = Value, color = Type), shape = 4) +
    labs(y = expression("OD"[450]), x = "%", fill = "") +
    scale_fill_identity(guide = "legend", labels = c("Specificity", "Sensitivity")) +
    theme_bw() + theme(aspect.ratio = 0.7) + theme(legend.position = "bottom") +
    facet_grid(cols = vars(protein)) +
    geom_text(data = df_label, aes(x = x, y = y), label = c(expression("f"["sp"]^"S"),
                                                            expression("f"["se"]^"S"),
                                                            expression("f"["sp"]^"N"),
                                                            expression("f"["se"]^"N")))
  ggsave(filename = here::here("outputs", "figs", "supp_sens_spec_od.pdf"),
    height = height, width = width)
}

plot_roc <- function(height, width) {
  datafit_roc <- get_datafit_control() %>%
    select(od, Value, Type, protein) %>%
    tidyr::spread(Type, Value) %>%
    mutate(FNR = (100 - Sensitivity) / 100, FPR = (100 - Specificity) / 100)
  roc_plot <- data.frame()
  for (i in seq(0, 2, 0.2)) {
    roc_plot <- rbind(roc_plot,
      filter(datafit_roc, protein == "N")[which.min(abs(filter(datafit_roc, protein == "N")$od - i)), ])
    roc_plot <- rbind(roc_plot,
      filter(datafit_roc, protein == "S")[which.min(abs(filter(datafit_roc, protein == "S")$od - i)), ])
  }
  roc_plot$od_label <-  seq(0, 2, 0.2) %>% purrr::map(~rep(.x, 2)) %>% unlist

  ggplot(data = datafit_roc, aes(x = FPR, y = Sensitivity / 100)) +
    geom_line(shape = 4) + geom_abline() +
    theme_bw() + theme(legend.position = "bottom") +
    geom_point(data = roc_plot, aes(y = Sensitivity / 100, x = FPR), color = palette(rainbow(22))) +
    geom_text(data = roc_plot, aes(y = Sensitivity / 100 - 0.024, x = FPR + 0.04, label = od_label),
      color = "red", size = 2.5) +
    facet_grid(cols = vars(protein)) + theme(aspect.ratio = 0.8)
  ggsave(filename = here::here("outputs", "figs", "roc.pdf"),
    height = height, width = width)
}

plot_roc_age <- function(height, width) {
  datafit_roc <- get_datafit_control() %>%
    select(od, Value, Type, protein) %>%
    tidyr::spread(Type, Value) %>%
    mutate(FNR = (100 - Sensitivity) / 100, FPR = (100 - Specificity) / 100)

  datafit_raw <- data_full$full %>%
    select(IgG_S_OD, IgG_N_OD, age_group, sero_pos)
  base_sens <- datafit_raw %>%
    filter(sero_pos == 1) %>%
    group_by(age_group) %>%
    summarise(freq = sum(sero_pos)) %>%
    pull(freq)
  base_spec <- datafit_raw %>%
    filter(sero_pos == 0) %>%
    group_by(age_group) %>%
    mutate(sero_pos = 1) %>%
    summarise(freq = sum(sero_pos)) %>%
    pull(freq)

  datafit_full <- data.frame()
  for (i in seq(0, 3, 0.01)) {
    odsum_sens <-  datafit_raw %>%
      filter(sero_pos == 1) %>%
      mutate(sero_pos_pred = (i < IgG_S_OD)) %>%
      group_by(age_group) %>%
      summarise(freq = sum(sero_pos_pred)) %>%
      mutate(sens = freq / base_sens, OD = i, protein = "S") %>%
      select(age_group, sens, protein, OD)
    odsum_spec <- datafit_raw %>%
      filter(sero_pos == 0) %>%
      mutate(sero_pos_pred = (i > IgG_S_OD)) %>%
      group_by(age_group) %>%
      summarise(freq = sum(sero_pos_pred)) %>%
      mutate(spec = freq / base_spec, OD = i, protein = "S") %>%
      select(age_group, spec, protein, OD)
    odsum_s <- merge(odsum_sens, odsum_spec, by = c("age_group", "OD", "protein"))

    odsum_sens <-  datafit_raw %>%
      filter(sero_pos == 1) %>%
      mutate(sero_pos_pred = (i < IgG_N_OD)) %>%
      group_by(age_group) %>%
      summarise(freq = sum(sero_pos_pred)) %>%
      mutate(sens = freq / base_sens, OD = i, protein = "N") %>%
      select(age_group, sens, protein, OD)
    odsum_spec <- datafit_raw %>%
      filter(sero_pos == 0) %>%
      mutate(sero_pos_pred = (i > IgG_N_OD)) %>%
      group_by(age_group) %>%
      summarise(freq = sum(sero_pos_pred)) %>%
      mutate(spec = freq / base_spec, OD = i, protein = "N") %>%
      select(age_group, spec, protein, OD)
    odsum_n <- merge(odsum_sens, odsum_spec, by = c("age_group", "OD", "protein"))
    datafit_full <- rbind(datafit_full, odsum_s, odsum_n)
  }
  datafit_full %<>% select(age_group, sens, spec, protein, OD) %>%
    unique %>%
    mutate(FPR = 1 - spec) %>%
    arrange(sens)

  lol <- data.frame()
  for (i in seq(0, 2, 0.2)) {
    lol1 <- datafit_full %>% filter(OD == i)
    lol <- rbind(lol, lol1)
  }
  ggplot(data = datafit_full, aes(x = FPR, y = sens)) +
    geom_line(size = 0.5, alpha = 0.8) +
    theme_bw() + theme(legend.position = "bottom") +
    labs(x = "FPR", y = "Sensitivity") +
    facet_grid(rows = vars(protein), cols = vars(age_group)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      aspect.ratio = 1.5) +
    geom_point(data = lol, aes(y = sens, x = FPR), color = (palette(rainbow(9)) %>%
      purrr::map(~rep(.x, 10)) %>%
      unlist)) +
    geom_text(data = lol, aes(y = sens - 0.02, x = FPR + 0.02, label = OD),
      color = (palette(rainbow(9)) %>% purrr::map(~rep(.x, 10)) %>% unlist),
      size = 2.5, position = position_jitter(width = 0.003, height = 0.001)) +
    xlim(0, 0.2) + ylim(0.25, 1)
  ggsave(filename = here::here("outputs", "figs", "roc_age.pdf"),
    height = height, width = width)
}