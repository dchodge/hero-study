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
      od_val  <- datafit_control_pro[which.min(abs(datafit_control_pro$Value - (i + 59))), ]$od
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

  datafitgp <- obj$datafit %>%
    group_by(!!arg_name, age_group, protein) %>%
    summarise(pos = sum(ind), n = n()) %>%
    mutate(!!arg_name := !!arg_name - 59)

  posti <- list()
  props <- list()
  for (pro in c("S", "N")) {
    datafitgp_pro <- datafitgp %>% filter(protein == pro)
    datalist <- list(
      N_gp = length(unique(eval(expr(`$`(datafitgp_pro, !!(arg_name)))))),
      x_gp = unique(eval(expr(`$`(datafitgp_pro, !!(arg_name))))),
      N = length(unique(eval(expr(`$`(datafitgp_pro, !!(arg_name)))))),
      D = length(unique(datafitgp_pro$age_group)),
      n = datafitgp_pro$n %>% matrix(, 5, byrow = TRUE),
      y = datafitgp_pro$pos %>% matrix(, 5, byrow = TRUE)
    )
    mod <- cmdstan_model(here::here("include", "multicontrol.stan"))
    fit <- mod$sample(data = datalist,
      chains = 4,
      parallel_chains = 4)

    posti[[pro]] <- fit$draws("p_post") %>%
      as_draws_df %>%
      spread_draws(p_post[s, a]) %>%
      ggdist::median_qi() %>%
      mutate(pro = pro, s = s + 59)
  }
  posti_all <- bind_rows(posti)

  eval(rlang::parse_expr(paste0(obj$fit_name, " <- posti_all")))
  outsting <- paste0("fit_mcmc_", obj$short_name, ".RData")
  saveexp <- substitute(save(x, file = y),
                        list(x =  obj$fit_name,
                             y =  rlang::expr(here::here("data", !!outsting))))
  eval(saveexp)

}

#' @export
plot <- function(obj, width, height, file_name) UseMethod("plot", obj)
#' @export
plot.test_metric <- function(obj, width, height, file_name) {
 # if (!exists(obj$fit_name)) {
    exec("data",  paste0("fit_mcmc_", obj$short_name))
 # }
  tots_age <- (obj$datafit %>% group_by(age_group) %>% summarise(n = n()) %>% pull(n))
  props <- tots_age / sum(tots_age)

  posti_all <- eval(str2lang(obj$fit_name))
  posti_all$a <- factor(posti_all$a, levels = 1:5, labels = c("<30", "30-39", "40-49", "50-59", "60+"))
  posti_all$pro <- factor(posti_all$pro, levels = c("S", "N"), labels = c("spike", "NCP"))

  posti_all_mean <- posti_all %>%
    group_by(s, pro) %>%
    summarise(p_post = weighted.mean(p_post, props),
      .lower = weighted.mean(.lower, props),
      .upper = weighted.mean(.upper, props),
    )

  p1 <- ggplot(data = posti_all_mean) +
    geom_abline(colour = "gray") + 
    ggdist::geom_lineribbon(aes(x = s, y = p_post * 100, ymin = .lower * 100, ymax = .upper * 100),
      fill = "black", size = 1, alpha = 0.6) +
    facet_grid(rows = vars(pro)) + theme_bw() + theme(aspect.ratio = 0.3) +
    labs(x = paste0(obj$name, " of control (%)"), y = paste0("Implied ", tolower(obj$name), " of HERO study (%)"), 
      fill = "Age group")

  p2 <- ggplot(data = posti_all) +
    geom_abline(colour = "gray") + 
    ggdist::geom_lineribbon(aes(x = s, y = p_post * 100, ymin = .lower * 100, ymax = .upper * 100,
      fill = a), size = 1, alpha = 0.6) +
    facet_grid(rows = vars(pro)) + theme_bw() + theme(aspect.ratio = 0.3, legend.position = "top") +
    labs(x = paste0(obj$name, " of control (%)"), 
      y =  paste0("Implied age-specific ", tolower(obj$name), "\n of HERO study (%)"), fill = "Age group")
  p1 / p2 + plot_annotation(tag_levels = 'A')
  ggsave(filename = here::here("outputs", "figs", paste0(file_name, ".pdf")), width = width, height = height)
}