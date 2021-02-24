plot_threshold_var <- function(thres_S, thres_N) {
  datafit_OD_val_i <- get_datafit_od_val(thres_S, thres_N)
  sens_S <- filter(datafit_OD_val_i, ab_protein == "S", Pos_Neg == "Pos") %>%
    select(sens_ind) %>%
    summarise(s = mean(sens_ind) * 100) %>%
    round(digits = 2) %>%
    as.numeric()
  spec_S <- filter(datafit_OD_val_i, ab_protein == "S", Pos_Neg == "Neg") %>%
    select(sens_ind) %>%
    summarise(s = mean(sens_ind) * 100) %>%
    round(digits = 2) %>%
    as.numeric()
  sens_N <- filter(datafit_OD_val_i, ab_protein == "N", Pos_Neg == "Pos") %>%
    select(sens_ind) %>%
    summarise(s = mean(sens_ind) * 100) %>%
    round(digits = 2) %>%
    as.numeric()
  spec_N <- filter(datafit_OD_val_i, ab_protein == "N", Pos_Neg == "Neg") %>%
    select(sens_ind) %>%
    summarise(s = mean(sens_ind) * 100) %>%
    round(digits = 2) %>%
    as.numeric()
  ggplot(data = datafit_OD_val_i) +
    geom_hline(aes(yintercept = od_thres), linetype = "dashed") +
    geom_jitter(aes(x = Pos_Neg, y = OD, color = sens_ind), shape = 21, alpha = 0.9, size = 2, width = 0.3) +
    theme_bw() +
    labs(x = "Known serostatus of individual", y = expression("OD"[450]), color = "Correctly identified?") +
    theme(legend.position = "top") +
    facet_grid(cols = vars(ab_protein)) +
    geom_text(aes(x, y, label = paste0("Sens: ", lab)),
      data = data.frame(x = 1, y = 4, lab = c(sens_S, sens_N),
        ab_protein = c("S", "N")), vjust = 1) +
    geom_text(aes(x, y, label = paste0("Spec: ", lab)),
      data = data.frame(x = 2, y = 4, lab = c(spec_S, spec_N),
        ab_protein = c("S", "N")), vjust = 1)
}