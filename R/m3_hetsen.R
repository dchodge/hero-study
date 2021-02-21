#### #### #### #### #### #### #### #### ####
#### AGE-SPECIFIC SENS + SPEC FROM CONTROL DATA ####
#### #### #### #### #### #### #### #### ####
# Class for the sens and spec
#' @export
test_metric <- function(name, vals, dataraw, func) {
  obj <- list()
  obj$name <- name
  obj$short_name <- rlang::parse_expr(tolower(substr(name, 1, 4)))
  obj$fit_name <- paste0("m3_", obj$short_name)
  obj$vals <- vals
  obj$func <- func
  obj$dataraw <- dataraw
  obj$datafit <- get_datafit_metric(obj)
  class(obj) <- "test_metric"
  obj
}

#' function to clean the data in object
get_datafit_metric <- function(obj) {
  dffit <- vector(mode = "list", length = 2 * length(obj$vals))
  j <- 1
  for (pro in c("S", "N")) {
    datafit_control_pro <- datafit_control %>% filter(protein == pro, Type == obj$name)
    for (i in 1:41) {
      od_val  <- datafit_control_pro[which.min(abs(datafit_control_pro$Value - (i + 59))), ]$OD
      dffit[[j]] <- data.frame(sample = seq_len(nrow(obj$dataraw)),
                               ind = obj$func(obj$dataraw[[paste0("IgG_", pro, "_OD")]], od_val) %>% as.numeric,
                               age_group = obj$dataraw$age_group,
                               temp = obj$vals[i],
                               protein = pro) %>%
                               rename(!!obj$short_name := temp)
      j <- j + 1
    }
  }
  data_reg <- bind_rows(dffit) %>%
    mutate(age_group = fct_recode(age_group,
      "1" = levels(obj$dataraw$age_group)[1],
      "2" = levels(obj$dataraw$age_group)[2], "3" = levels(obj$dataraw$age_group)[3],
      "4" = levels(obj$dataraw$age_group)[4], "5" = levels(obj$dataraw$age_group)[5]))
  data_reg
}

#' @export
fit <- function(obj) UseMethod("fit", obj)
#' @export
fit.test_metric <- function(obj) {
  arg_name <- obj$short_name
  arg_val <- as.character(obj$vals)
  form <- paste0("ind ~ (1 | ", arg_name, "/age_group)")
  fitted <- c("S", "N") %>%
    map(~stan_glmer(form, data = dplyr::filter(obj$datafit, protein == .x), family = binomial(link = "logit")))
  eval(rlang::parse_expr(paste0(obj$fit_name, " <- fitted")))
  outsting <- paste0("fit_mcmc_", obj$short_name, ".RData")
  saveexp <- substitute(save(x, file = y),
                        list(x =  obj$fit_name,
                             y =  rlang::expr(here::here("data", !!outsting))))
  eval(saveexp)
}

#' @export
plot <- function(obj) UseMethod("plot", obj)
#' @export
plot.test_metric <- function(obj) {
  if (!exists(obj$fit_name)) {
    exec("data",  paste0("fit_mcmc_", obj$short_name))
  }
  arg_name <- obj$short_name
  arg_val <- as.character(obj$vals)
  df_pred <- exec("expand.grid", age_group = 1:5, !!arg_name := arg_val)
  arg1 <- rlang::parse_expr(paste0(obj$fit_name, "[[1]]"))
  arg2 <- rlang::parse_expr(paste0(obj$fit_name, "[[2]]"))
  df_pred_both <- rbind(eval(rlang::call2("add_fitted_draws", df_pred, arg1, n = 100)) %>%
    mutate(protein = "S", measure = !!obj$name),
                        eval(rlang::call2("add_fitted_draws", df_pred, arg2, n = 100)) %>%
                        mutate(protein = "N", measure = !!obj$name))
  df_pred_both$age_group <- factor(df_pred_both$age_group, levels = 1:5,
    labels = c("<30", "30-39", "40-49", "50-59", "60+"))
  p1 <- df_pred_both %>%
    ggplot(aes(x = as.numeric(!!obj$short_name) + 59)) +
    geom_line(aes(y = .value * 100, group = .draw), alpha = 0.05, color = "#08519C")  +
    stat_summary(aes(y = .value * 100), fun = median,  geom = "line", color = "red", size = 0.75) +
    facet_grid(cols = vars(age_group), rows = vars(protein)) +
    labs(x = paste0(obj$name, " of control"),
         y =  paste0("Implied ", obj$name)) +
    theme_bw()
  p2 <- df_pred_both %>%
    ggplot(aes(x = as.numeric(!!obj$short_name) + 59, color = age_group)) +
    stat_summary(aes(y = .value * 100, group = age_group), fun = median,  geom = "line",  size = 0.75) +
    theme_bw() + theme(legend.position = "top") + labs(x = "Sensitivity of control", y = "Implied sensitivity") +
    facet_grid(rows = vars(protein)) +
    labs(x = paste0(obj$name, " of control"),
         y =  paste0("Implied ", obj$name),
         color = "Age group (yrs)")
  ggpubr::ggarrange(p1, p2, nrow = 2)
}