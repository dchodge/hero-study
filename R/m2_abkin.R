#### #### #### #### #### #### #### #### ####
#### Predict titre levels at first bleed ####
#### #### #### #### #### #### #### #### ####
# Class to assist seroprevalence fitting
#' @export
init_model_ab_kin <- function(name, stan_file, datafit, marginal_vars, y_axis) {
  obj <- list()
  obj$name <- name
  obj$fit_name <- paste0("m2_", name)
  obj$stan_file <- here::here("include", paste0(stan_file, ".stan"))
  obj$stan_file_out <- here::here(paste0("include/", stan_file, "/"))
  obj$datafit <- datafit
  obj$marginal <- marginal_vars
  obj$y_axis <- y_axis
  obj$measure <- c("IgG_S", "IgG_N", "IgA_S_Serum", "IgA_N_Serum")
  obj$model_list <- list()
  class(obj) <- "Ab_model"
  obj
}

#' @export
fit <- function(Ab_model) UseMethod("fit", Ab_model)
#' @export
fit.Ab_model <- function(obj) {
  model <- cmdstan_model(obj$stan_file, compile = FALSE)
  m_compile <- model$compile(dir = obj$stan_file_out, force_recompile = TRUE)

  data_stan_format <- obj$datafit %>%
    select(age_group, gender, ethnic, symp_pos) %>%
    tidybayes::compose_data()

  obj$model_list <- list()
  for (i in 1:4) {
    data_list <- list(
      n = data_stan_format$n,
      num_age = data_stan_format$n_age_group,
      num_eth = data_stan_format$n_ethnic,
      num_gen = data_stan_format$n_gender,
      num_symp = data_stan_format$n_symp_pos,
      age_group = data_stan_format$age_group,
      eth = data_stan_format$ethnic,
      gen = data_stan_format$gender,
      symp = data_stan_format$symp_pos,
      p_eth = as.numeric((obj$datafit$ethnic %>% table) / nrow(obj$datafit)),
      p_gen = as.numeric((obj$datafit$gender %>% table) / nrow(obj$datafit)),
      p_age = as.numeric((obj$datafit$age_group %>% table) / nrow(obj$datafit)),
      p_symp = as.numeric((obj$datafit$symp_pos %>% table) / nrow(obj$datafit)),
      ab_measure = obj$datafit[[obj$measure[i]]],
      Time = obj$datafit$time
    )
    fit <- m_compile$sample(
      data = data_list,    # named list of data
      iter_warmup = 1000,          # number of warmup iterations per chain
      iter_sampling = 1000,            # total number of iterations per chain
      chains = 4,            # number of cores (could use one per chain)
      parallel_chains = 4,
      adapt_delta = 0.999,
      max_treedepth = 20
    )
    obj$model_list[[i]] <- fit$draws()
  }
  # saving
  eval(rlang::parse_expr(paste0(obj$fit_name, " <- obj$model_list")))
  outsting <- paste0("fit_mcmc_", obj$name, ".RData")
  saveexp <- substitute(save(x, file = y),
                        list(x =  obj$fit_name,
                             y =  rlang::expr(here::here("data", !!outsting))))
  eval(saveexp)
}

#' @export
clean_samples.Ab_model <- function(Ab_model) UseMethod("clean_samples", Ab_model)
#' @export
clean_samples.Ab_model <- function(obj) {
  if (!exists(obj$fit_name)) {
    dataexpr <- expr(data(source))
    dataexpr[[2]] <- paste0("fit_mcmc_", obj$name)
    eval(dataexpr)
  }
  dataexpr <- expr(data(source))
  dataexpr[[2]] <- paste0("fit_mcmc_", obj$name)
  eval(dataexpr)
  cov_names <- c("age_group", "ethnic", "gender", "symp_pos")
  variables <- cov_names %>% purrr::map(~levels(get_datafit_start(FALSE)[[.x]]))

  df_plot_lol <- data.frame()
  for (m in 1:4) {
    # extracts the data from obj$fit_name"[[m]]" and saves in posterior
    z <- expr(posteriors <- fit_name)
    z[[3]] <- rlang::parse_expr(paste0(obj$fit_name, "[[m]]"))
    eval(z)
    for (i in 1:4) {
      x <- str2lang(paste0(obj$marginal[i], "[j]"))
      df_temp <- posteriors %>%
        as_draws_df %>%
        tidybayes::spread_draws(!!x) %>%
        ggdist::median_qi() %>%
        as.data.frame %>%
        select(`50%` := !!obj$marginal[i], `2.5%` = .lower, `97.5%` = .upper)
      df_temp$cov_names <- cov_names[i]
      df_temp$cov <- variables[[i]]
      df_temp$measure <- obj$measure[m]
      df_plot_lol <- rbind(df_plot_lol, df_temp)
    }
  }
  df_plot <- df_plot_lol
  df_plot
}

clean_data_for_plt <- function() {
  df_long_temp_start_long <- get_datafit_start() %>%
    select(!c(Sample_ID, time, date_of_sample, sample_no, days_since_inf)) %>%
    tidyr::gather(key  = measure, value, -ethnic, -gender, -age_group, -symp_pos) %>%
    tidyr::gather(key  = cov_names, cov, -measure, -value) %>%
    group_by(measure, cov, cov_names) %>%
      na.omit
  df_long_temp_start_long <- df_long_temp_start_long %>%
    dplyr::summarise(value = mean(value)) %>%
    dplyr::arrange(order(match(measure, df_long_temp_start_long$measure)))
  df_long_temp_change_long <- get_datafit_change() %>%
    select(!c(time)) %>%
    tidyr::gather(key  = measure, value, -ethnic, -gender, -age_group, -symp_pos) %>%
    tidyr::gather(key  = cov_names, cov, -measure, -value) %>%
    group_by(measure, cov, cov_names) %>%
    na.omit
  df_long_temp_change_long <- df_long_temp_change_long %>%
    dplyr::summarise(value = mean(value)) %>%
    dplyr::arrange(order(match(measure, df_long_temp_change_long$measure)))
  plt_data <- list(start = df_long_temp_start_long, change = df_long_temp_change_long)
  plt_data
}

#' @export
plot.Ab_model <- function(Ab_model) UseMethod("plot", Ab_model)
#' @export
plot.Ab_model <- function(obj) {
  measures <- c("IgG_S", "IgG_N", "IgA_S_Serum", "IgA_N_Serum")
  measures_name <- c("spike-IgG", "NCP-IgG", "spike-IgA", "NCP-IgA")

  cov_names_1 <- c("age_group", "ethnic", "gender", "symp_pos")
  cov_names_lab <- c("Age group", "Ethnicity", "Gender", "Disease severity")
  dataplt <- clean_samples(obj)
  dataplt <- dataplt %>%
    mutate(measure = factor(measure, levels = measures, labels = measures_name)) %>%
    mutate(cov_names = factor(cov_names, levels = cov_names_1, labels = cov_names_lab))
  get_plot_data <- clean_data_for_plt()[[obj$name]] %>%
    mutate(measure = factor(measure, levels = measures, labels = measures)) %>%
    mutate(cov_names = factor(cov_names, levels = cov_names_1, labels = cov_names_lab))
  ggplot(data = dataplt) +
  #  geom_point(position = position_dodge(width = 0.5),
  #    data = get_plot_data, aes(y = value, x = cov,
  #    color = factor(measure, levels = measures, labels = measures)), shape = 4, size = 0.8) +
    geom_pointrange(aes(ymin = `2.5%`, y = `50%`, ymax = `97.5%`, x = cov, fill = measure, color = measure),
      position = position_dodge(width = 0.5), alpha = 0.5) +
    facet_grid(cols = vars(cov_names), scale = "free_x") +
#    stat_summary(aes(y = `50%`, x = cov, group = measure, color = measure), fun = mean,
#      geom = "line",  position = position_dodge(width = 0.5)) +
   ## stat_summary(aes(y = `50%`, x = cov, group = measure, color = measure), fun = mean,
  #    geom = "point", shape = 5,  position = position_dodge(width = 0.5)) +
    theme_bw() + labs(y = obj$y_axis, x = "", colour = "") +
    theme(axis.title.x = element_blank(), axis.title.y = element_text(size = 10),
      axis.text.x = element_blank()) +
    scale_color_manual(values = c("#b31818", "#e95d5d", "#00adcc", "#1eddff")) +
    scale_fill_manual(values = c("#b31818", "#e95d5d", "#00adcc", "#1eddff"), guide = FALSE)
}


clean_samples_revert <- function(obj_start, obj_change) {
  dataexpr <- expr(data(source))
  dataexpr[[2]] <- paste0("fit_mcmc_", obj_start$name)
  eval(dataexpr)

  dataexpr <- expr(data(source))
  dataexpr[[2]] <- paste0("fit_mcmc_", obj_change$name)
  eval(dataexpr)
  
  cov_names <- c("age_group", "ethnic", "gender", "symp_pos")
  variables <- cov_names %>% purrr::map(~levels(get_datafit_start()[[.x]]))
  df_plot_lol <- data.frame()
  for (m in 1:4) {
    revert_val <- min(obj_start$datafit[[obj_start$measure[m]]])
    z1 <- expr(post_s <- fit_name)
    z1[[3]] <- rlang::parse_expr(paste0(obj_start$fit_name, "[[m]]"))
    eval(z1)
    z2 <- expr(post_c <- fit_name)
    z2[[3]] <- rlang::parse_expr(paste0(obj_change$fit_name, "[[m]]"))
    eval(z2)

    for (i in 1:4) {
      post_marg_s <- post_s %>% as_draws_df %>% select(starts_with(obj_start$marginal[i]))
      post_marg_c <- post_c %>% as_draws_df %>% select(starts_with(obj_change$marginal[i]))
      time_revert <- (revert_val - post_marg_s) * 28 / post_marg_c / 7
      df_temp <- time_revert %>% apply(2, quantile, c(0.025, 0.5, 0.975)) %>% t %>% as.data.frame
      df_temp$cov_names <- cov_names[i]
      df_temp$cov <- variables[[i]]
      df_temp$measure <- obj_start$measure[m]
      df_plot_lol <- rbind(df_plot_lol, df_temp)
    }
  }
  df_plot_revert <- df_plot_lol
  df_plot_revert
}

plot_revert  <- function(obj_start, obj_change) {
  dataplt_revert <- clean_samples_revert(obj_start, obj_change)
  measures <- c("IgG_S", "IgG_N", "IgA_S_Serum", "IgA_N_Serum")
  cov_names_1 <- c("age_group", "ethnic", "gender", "symp_pos")
  cov_names_lab <- c("Age group", "Ethnicity", "Gender", "Disease severity")
  dataplt_revert$measure <- factor(dataplt_revert$measure, levels = measures, labels = measures)
  dataplt_revert$cov_names <- factor(dataplt_revert$cov_names, levels = cov_names_1, labels = cov_names_lab)

  ggplot(data = dataplt_revert) +
    geom_pointrange(aes(ymin = `2.5%`, y = `50%`, ymax = `97.5%`, x = cov, color = measure),
      position = position_dodge(width = 0.5), alpha = 0.5) +
    facet_grid(cols = vars(cov_names), scale = "free_x") +
  #  stat_summary(aes(y = `50%`, x = cov, group = measure, color = measure), fun = mean,
  #    geom = "line",  position = position_dodge(width = 0.5)) +
  #  stat_summary(aes(y = `50%`, x = cov, group = measure, color = measure), fun = mean,
 #     geom = "point", shape = 5,  position = position_dodge(width = 0.5)) +
    theme_bw() + labs(y = "Week PSO until seroreversion", x = "", color = "") +
    theme(axis.title.y = element_text(size = 10), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_color_manual(values = c("#b31818", "#e95d5d", "#00adcc", "#1eddff")) +
      coord_cartesian(ylim = c(0, 200))
}

plot_figure_2 <- function(obj_start, obj_change, width_i =  10, height_i = 10) {
  p1 <- plot(obj_start)
  p2 <- plot(obj_change)
  p3 <- plot_revert(obj_start, obj_change)
  ggpubr::ggarrange(p1, p2, p3, nrow = 3, common.legend = TRUE)
  ggsave(filename = here::here("outputs", "figs", "fig2.pdf"), width = width_i, height = height_i)
}
