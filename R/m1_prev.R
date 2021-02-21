# Class to assist seroprevalence fitting
#' @export
init_model_prev <-  function(name, response, datafit, stan_file) {
  obj <- list()
  obj$cov <- c()
  obj$response <- response
  obj$stan_file <- here::here("include", paste0(stan_file, ".stan"))
  obj$stan_file_out <- here::here(paste0("include/", stan_file, "/"))
  obj$name <- name
  obj$datafit <- datafit
  obj$fit_name <- paste0("m1_", name)
  obj$model_col_name <- "unnamed"
  obj$output_summary <- TRUE
  obj$covar_labels <- list(
    c19zones = c("1", "2", "3", "4", "5", "6"),
    job_grouped = c(datafit$job_grouped %>% levels),
    loc_coded = c(datafit$loc_coded %>% levels),
    gender = c(datafit$gender %>% levels),
    ethnic = c(datafit$ethnic %>% levels),
    age_group  = datafit_prev$age_group %>% levels
  )
  obj$pri_covar_l <- list("c19zones", "job_grouped", "loc_coded")
  obj$model_num <- list(c19zones = "Model A", job_grouped = "Model B", loc_coded = "Model C")
  class(obj) <- "Prev_model"
  obj
}

# Fit the seroprevalence model for a given class instance
sampling_helper <- function(obj, stan_model, file_name) {

  data_stan <- tidybayes::compose_data(obj$datafit %>% select(c(obj$cov, obj$response)))
  data_convert <- data_stan[c(obj$cov, obj$response)] %>%
    as.data.frame %>%
    dplyr::group_by_at(vars(one_of(obj$cov))) %>%
    dplyr::summarise(response = sum(!!sym(obj$response)), n = n())
  cov_size <- obj$cov %>% map_dbl(~data_stan[[paste0("n_", .x)]]) %>% prepend(0)
  indicies <- obj$cov %>% map(~1:data_stan[[paste0("n_", .x)]]) %>% expand.grid
  max.len <- max(cov_size)
  prop <- list()
  for (covar in obj$cov) {
    tabs <- table(data_convert[[covar]]) / length(data_convert[[covar]])
    prop[[covar]] <- c(tabs, rep(Inf, max.len - length(tabs)))
  }
  prop <- as.data.frame(prop)

  datalist <- list(
    n = data_convert %>% nrow,
    n_cov = length(obj$cov),
    y_se = 179,
    n_se = 180,
    y_sp = 646,
    n_sp = 650,
    response = data_convert$response,
    tot = data_convert$n,
    covariates = data_convert[obj$cov],
    cov_size = cov_size,
    cov_size_tot = sum(cov_size),
    indicies = indicies,
    prop = prop
  )
  m1 <- sampling(stan_model, data = datalist,
                 chains = 4, cores = 4,
                 control = list(adapt_delta = 0.999, max_treedepth = 20),
                 sample_file = paste0(obj$stan_file_out, file_name),
                 diagnostic_file = paste0(obj$stan_file_out, file_name, "_diag")
                 )
  m1
}

#' @export
fit <- function(Prev_model) UseMethod("fit", Prev_model)
#' @export
fit.Prev_model <- function(obj) {
  stan_sero <- stan_model(file = obj$stan_file)
  obj$fit <- list()
  for (pri_covar in obj$pri_covar) {
    obj$cov <-  c(pri_covar)
    obj$fit[[pri_covar]][["unadjusted"]] <- sampling_helper(obj, stan_sero, paste0(pri_covar, "_unadjusted"))
    obj$cov <-  c(pri_covar, "gender", "ethnic", "age_group")
    obj$fit[[pri_covar]][["adjusted"]] <- sampling_helper(obj, stan_sero, paste0(pri_covar, "_adjusted"))
  }
  eval(rlang::parse_expr(paste0(obj$fit_name, " <- obj$fit"))) # save obj$fit to fitname
  outsting <- paste0("fit_mcmc_", obj$name, ".RData") # name of file
  saveexp <- substitute(save(x, file = y),
                        list(x =  rlang::parse_expr(obj$fit_name),
                             y = rlang::expr(here::here("data", !!outsting)))) # create save expression
  eval(saveexp) # execute save
}

#' @export
clean_samples <- function(Prev_model) UseMethod("clean_samples", Prev_model)
#' @export
clean_samples.Prev_model <- function(obj) {
  if (!exists(obj$fit_name)) {
    dataexpr <- expr(data(source))
    dataexpr[[2]] <- paste0("fit_mcmc_", obj$name)
    eval(dataexpr)
  }
  df_post <- data.frame()
  df_post_raw <- data.frame()
  for (pri_covar in obj$pri_covar_l) {
    x_unad_expr <- rlang::parse_expr(paste0(obj$fit_name, "[[pri_covar]]$unadjusted"))
    x_ad_expr <- rlang::parse_expr(paste0(obj$fit_name, "[[pri_covar]]$adjusted"))
    df_post_un <- tidybayes::spread_draws(eval(x_unad_expr), marginal[group, cov]) %>%
      mutate(cov = case_when(cov == 1~pri_covar), type = "unadjusted", model =  obj$model_num[[pri_covar]])
    df_post_ad <- tidybayes::spread_draws(eval(x_ad_expr), marginal[group, cov]) %>%
      mutate(cov = case_when(cov == 1~pri_covar, cov == 2~"gender", cov == 3~"ethnic", cov == 4~"age_group"),
             type = "adjusted", model =  obj$model_num[[pri_covar]])
    df_post_raw <- rbind(df_post_un, df_post_ad, df_post_raw)
  }
  df_post_raw %<>% group_by(cov, group, type, model) %>%
    summarise(mutate_quantile(marginal)) %>%
    tidyr::spread(probs, value)
  df_data_raw <- data.frame()
  obj$datafit$c19zones <- as.character(as.numeric(obj$datafit$c19zones))
  for (cov in names(obj$covar_labels)) {
    df_post_t <- df_post_raw %>%
      filter(cov == cov) %>%
      mutate(group = factor(group, levels = seq_len(obj$covar_labels[[cov]]), labels = obj$covar_labels[[cov]]))
    df_post <- rbind(df_post, df_post_t)
    temp <- obj$datafit %>%
      dplyr::group_by_at(cov) %>%
      dplyr::summarise(prev = mean(!!sym(obj$response))) %>%
      setNames(c("cat", "value")) %>%
      mutate(cov = cov)
    df_data_raw <- rbind(df_data_raw, temp)
  }
  obs_data <- rbind(
    df_data_raw %>%
      dplyr::filter(cov %in% c("c19zones", "gender", "ethnic", "age_group")) %>%
      mutate(model = "Model A"),
    df_data_raw %>%
      dplyr::filter(cov %in% c("job_grouped", "gender", "ethnic", "age_group")) %>%
      mutate(model = "Model B"),
    df_data_raw %>%
      dplyr::filter(cov %in% c("loc_coded", "gender", "ethnic", "age_group")) %>%
      mutate(model = "Model C")
  )
  # Add factors to covariates
  obs_data$cov <- factor(obs_data$cov, levels = colnames(datafit_prev)[1:6], labels = colnames(datafit_prev)[1:6])
  df_post$cov <- factor(df_post$cov, levels = colnames(datafit_prev)[1:6], labels = colnames(datafit_prev)[1:6])
  fig_data <- list(fit = df_post, data = obs_data)
  fig_data
}


#' @export
plot <- function(Prev_model, x_axis_str, file_name, height, width) UseMethod("plot", Prev_model)
#' @export
plot.Prev_model <- function(obj, x_axis_str, file_name, height, width) {
  plt_data <- clean_samples(obj)
  plt_data$fit$cov <- factor(plt_data$fit$cov, levels = levels(plt_data$fit$cov),
                             labels = c("Zones", "Job", "Job location", "Gen.",  "Age", "Ethnicity"))
  plt_data$data$cov <- factor(plt_data$data$cov, levels = levels(plt_data$data$cov),
                             labels = c("Zones", "Job", "Job location", "Gen.", "Age", "Ethnicity"))
  p1 <- ggplot(data = plt_data$fit) +
    geom_point(data = plt_data$data, aes(x = value, y = cat), color = "black", shape = 8, alpha = 1, size = 3) +
    geom_linerange(aes(xmin = `0.025`, xmax = `0.975`, y = group, colour = type),
      position = position_dodge(width = 0.9), size = 0.80) +
    geom_point(aes(x = `0.5`, y = group, fill = type), position = position_dodge(width = 0.9),
      shape = 21, size = 2, alpha = 0.8, show.legend = F) +
    facet_grid(rows = vars(cov), cols = vars(model), scales = "free_y", space = "free_y") +
    labs(x = x_axis_str, y = "Covariates", color = "Regression type", fill = "") +
    theme_bw() + theme(legend.position = "top")
  ggsave(filename = here::here("outputs", "figs", paste0(file_name, ".pdf")), height = height, width = width)
  p1
}