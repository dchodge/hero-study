# Generate datasets for models
get_datafit_prev <- function() {
  df_prev <- data_full$full %>%
    dplyr::filter(sample_no == 1) %>%
    dplyr::select(c19zones, job_grouped, loc_coded, gender, age_group, ethnic, sero_pos) %>%
    na.omit %>%
    tibble
  df_prev
}

get_datafit_asymp <- function() {
  df_prev <- data_full$full %>%
    dplyr::filter(sample_no == 1) %>%
    dplyr::select(c19zones, job_grouped, loc_coded, gender, age_group, ethnic, sero_pos, symp_pos) %>%
    na.omit %>%
    tibble
  df_prev_asymp <- df_prev %>%
  dplyr::filter(sero_pos == 1) %>%
    mutate(job_grouped = forcats::fct_drop(job_grouped), asymp_pos = abs(symp_pos - 1)) %>%
    tibble

  df_prev_asymp
}

get_datafit_start <- function(numerise = TRUE) {
  df_ab_long_d <- data_full$sero_pos %>%
    dplyr::select(Sample_ID,
                  IgG_S = IgG_S_unitl2_who,
                  IgG_N = IgG_N_unitl2_who,
                  IgA_S_Serum = IgA_Serum_N_unitl2,
                  IgA_N_Serum = IgA_Serum_S_unitl2,
                  age_group,  gender, ethnic, symp_pos, date_of_sample, sample_no, days_since_inf) %>%
    dplyr::filter(sample_no == 1) %>%
      mutate(time = 1, symp_pos = factor(symp_pos, levels = c(0, 1), labels = c("Asymp", "Symp")))
  df_ab_long_d
}

get_datafit_change <- function(numerise = TRUE) {
  df_ab_long_d <- data_full$sero_pos %>%
    dplyr::select(Sample_ID,
                  IgG_S = IgG_S_unitl2_who,
                  IgG_N = IgG_N_unitl2_who,
                  IgA_S_Serum = IgA_Serum_N_unitl2,
                  IgA_N_Serum = IgA_Serum_S_unitl2,
                  age_group,  gender, ethnic, symp_pos, date_of_sample, sample_no, days_since_inf) %>%
    mutate(symp_pos = factor(symp_pos, levels = c(0, 1), labels = c("Asymp", "Symp")))
  df_ab_long_d$time <- (df_ab_long_d$date_of_sample) -
    (purrr::map(df_ab_long_d$date_of_sample[seq(1, nrow(df_ab_long_d), 2)], ~rep(.x, 2)) %>% unlist)
  df_long_temp_start <- df_ab_long_d %>%
    dplyr::filter(sample_no == 1) %>%
    dplyr::select(IgG_S, IgG_N, IgA_S_Serum, IgA_N_Serum)
  df_long_temp_end <- df_ab_long_d %>%
    dplyr::filter(sample_no == 2) %>%
    dplyr::select(IgG_S, IgG_N, IgA_S_Serum, IgA_N_Serum)
  df_long_temp_change_num <- df_long_temp_end - df_long_temp_start

  df_long_temp_change <- dplyr::bind_cols(df_long_temp_change_num,
    dplyr::select(dplyr::filter(df_ab_long_d, sample_no == 2),
      time, age_group, ethnic, gender, symp_pos))
  df_long_temp_change
}

# Sort out level 3

get_datafit_seropos <- function() {
  hero_data <- data_full$full %>%
    dplyr::filter(sample_no == 1) %>%
    dplyr::select(Sample_ID, age, age_group, gender, ethnic, IgG_S_OD, IgG_N_OD, IgG_S_pos, IgG_N_pos, sero_pos) %>%
    filter(sero_pos == 1)
  hero_data
}

get_datafit_seroneg <- function() {
  hero_data <- data_full$full %>%
    dplyr::filter(sample_no == 1) %>%
    dplyr::select(Sample_ID, age, age_group, gender, ethnic, IgG_S_OD, IgG_N_OD, IgG_S_pos, IgG_N_pos, sero_pos) %>%
    filter(sero_pos == 0)
  hero_data
}

get_datafit_control <- function() {
  elisa_data_s_raw <- read.csv(system.file("extdata", "Sheffield_ELISA_validation_data_S.csv", package = "serobayes"))
  elisa_data_n_raw <- read.csv(system.file("extdata", "Sheffield_ELISA_validation_data_N.csv", package = "serobayes"))
  elisa_data_s <- clean_data_od_val_curve_protein(elisa_data_s_raw, "S")
  elisa_data_n <- clean_data_od_val_curve_protein(elisa_data_n_raw, "N")
  elisa_data <- rbind(elisa_data_s, elisa_data_n)
}

get_datafit_time <- function(numerise = TRUE) {
  df_ab_long_d <- data_full$symp %>%
    dplyr::select(Sample_ID,
                  IgG_S = IgG_S_unitl2_who,
                  IgG_N = IgG_N_unitl2_who,
                  IgA_S_Serum = IgA_Serum_N_unitl2,
                  IgA_N_Serum = IgA_Serum_S_unitl2,
                  age_group,  gender, ethnic, symp_pos, date_of_sample, sample_no, days_since_inf) %>%
      mutate(time = 1, symp_pos = factor(symp_pos, levels = c(0, 1), labels = c("Asymp", "Symp")))
  df_ab_long_d
}

clean_data_od_val_curve_protein <- function(elisa_data, protein_str) {
  od <- elisa_data$OD %>% purrr::map_dbl(~substr(.x, 3, nchar(.x)) %>% as.numeric)
  sens <- elisa_data$Sensitivity
  spec <- elisa_data$Specificity
  sens_ci_split <- unlist(strsplit(elisa_data$X95..CI, " "))
  sens_ci_lower <- sens_ci_split[seq(1, length(sens_ci_split), 3)] %>%
    purrr::map_dbl(~substr(.x, 1, nchar(.x) - 1) %>% as.numeric)
  sens_ci_upper <- sens_ci_split[seq(3, length(sens_ci_split), 3)] %>%
    purrr::map_dbl(~substr(.x, 1, nchar(.x) - 1) %>% as.numeric)
  spec_ci_split <- unlist(strsplit(elisa_data$X95..CI.1, " "))
  spec_ci_lower <- spec_ci_split[seq(1, length(spec_ci_split), 3)] %>%
    purrr::map_dbl(~substr(.x, 1, nchar(.x) - 1) %>% as.numeric)
  spec_ci_upper <- spec_ci_split[seq(3, length(spec_ci_split), 3)] %>%
    purrr::map_dbl(~substr(.x, 1, nchar(.x) - 1) %>% as.numeric)
  df_elisa <- rbind(
    data.frame(od = od,
      Value = sens, CI_lower = sens_ci_lower, CI_upper = sens_ci_upper, Type = "Sensitivity"),
    data.frame(od = od,
      Value = spec, CI_lower = spec_ci_lower, CI_upper = spec_ci_upper, Type = "Specificity"))
  df_elisa$protein <- protein_str
  df_elisa_trim <- df_elisa[df_elisa$od > 0, ]
  df_elisa_trim
}


get_datafit_od_val <- function(thres_s, thres_n) {
  val_s_data <- read.csv(system.file("extdata", "HERO_Validation_OD_21Dec2020_S.csv", package = "serobayes"))
  val_s_data$ab_protein <- "S"
  val_s_data$od_thres <- thres_s
  t1 <- val_s_data %>% dplyr::filter(Pos_Neg == "Pos") %>% dplyr::mutate(sens_ind = (OD > thres_s))
  t2 <- val_s_data %>% dplyr::filter(Pos_Neg == "Neg") %>% dplyr::mutate(sens_ind = (OD < thres_s))
  val_s_data <- rbind(t1, t2)

  val_n_data <- read.csv(system.file("extdata", "HERO_Validation_OD_21Dec2020_N.csv", package = "serobayes"))
  val_n_data$ab_protein <- "N"
  val_n_data$od_thres <- thres_n
  t1 <- val_n_data %>% dplyr::filter(Pos_Neg == "Pos") %>% dplyr::mutate(sens_ind = (OD > thres_n))
  t2 <- val_n_data %>% dplyr::filter(Pos_Neg == "Neg") %>% dplyr::mutate(sens_ind = (OD < thres_n))
  val_n_data <- rbind(t1, t2)
  val_data <- rbind(val_s_data, val_n_data)
  val_data$Pos_Neg <- factor(val_data$Pos_Neg, levels = c("Pos", "Neg"))
  val_data$ab_protein <- factor(val_data$ab_protein, levels = c("S", "N"))
  val_data
}

get_data_commercial <- function() {
  df_commerical_test_sens <- data.frame(
    manufactuer = c("Abbott", "Diasorin", "Vidas", "Roche", "Elisa", "Beckman", "Siemens"),
    protein_measure = c("N", "S", "S", "N", "S", "S", "S"),
    Val = c(84.7, 82.4, 89.3, 89.0, 89.4, 81.5, 85.9),
    CI_lower = c(81.5, 79.0, 85.5,  85.9, 85.4, 74.6, 79.4),
    CI_upper = c(87.4,  85.3, 92.1, 91.4, 92.3, 86.6, 90.4),
    Type = rep("Sensitivity", 7)
  )
  df_commerical_test_spec <- data.frame(
    manufactuer = c("Abbott", "Diasorin", "Vidas", "Roche", "Elisa", "Beckman", "Siemens"),
    protein_measure = c("N", "S", "S", "N", "S", "S", "S"),
    Val = c(99.5, 98.47, 98.9, 100, 97.7, 100, 99.8),
    CI_lower = c(99.2, 98.2, 98.2, 99.8, 96.7, 98.8, 98.7),
    CI_upper = c(99.74, 99.1, 99.4, 100, 98.4, 100, 99.94),
    Type = rep("Specificity", 7)
  )
  df_commerical_test <- rbind(df_commerical_test_sens, df_commerical_test_spec)
  df_commerical_test$type <- "observed"
  df_commerical_test
}