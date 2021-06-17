#' @import tidyr
#' @import dplyr
#' @import forcats
#' @import purrr
#' @import ggplot2
#' @import patchwork
#' @import posterior
#' @importFrom tidybayes add_fitted_draws compose_data spread_draws
#' @importFrom magrittr %>% %<>%
#' @importFrom rstan stan_model sampling extract
#' @importFrom rstanarm stan_glmer
NULL

# Pipelines for relabeling, cleaning, removing etc.
pipeline_rename <- function(data) {
  data %<>% rename(
         ethnic = EthnGroup,
         date_of_sample = date,
         IgG_S_OD = IgG_spike_OD,
         IgG_N_OD = IgG_NCP_OD,
         IgG_S_pos = IgG_spike_designation,
         IgG_N_pos = IgG_NCP_designation,
         IgG_S_titre = IgG_spike_endpoint,
         IgG_N_titre = IgG_NCP_endpoint,
         IgG_S_unit_hero = IgG_spike_antibody_units,
         IgG_N_unit_hero = IgG_NCP_antibody_units,
         sero_pos = serostatus,
         IgA_Serum_S_OD = Serum_IgA_spike_OD,
         IgA_Serum_S_unit = Serum_IgA_spike_antibody_units,
         IgA_Serum_N_OD = Serum_IgA_NCP_OD,
         IgA_Serum_N_unit = Serum_IgA_NCP_antibody_units
  )
  data
}


pipeline_relabel_cat <- function(data) {
  data_trans <- data %>%
    mutate(PCR_confirmed =
      dplyr::case_when(PCR_confirmed == ""~0, PCR_confirmed == "y"~1)) %>%
    mutate(IgG_S_pos =
      dplyr::case_when(IgG_S_pos == "Neg"~0, IgG_S_pos == "Pos"~1)) %>%
    mutate(IgG_N_pos =
      dplyr::case_when(IgG_N_pos == "Neg"~0, IgG_N_pos == "Pos"~1)) %>%
    mutate(sero_pos =
      dplyr::case_when(sero_pos == "Neg"~0, sero_pos == "Pos"~1)) %>%
    mutate(job_grouped =
      dplyr::case_when(job_grouped == "Nurse"~"nurse",
        job_grouped == "HCA"~"hca",
        job_grouped == "Pharmacy"~"pharmacy",
        job_grouped == "0"~"NA",
        TRUE ~ as.character(job_grouped))) %>%
    mutate(loc_coded = dplyr::case_when(loc_coded == "0"~"NA", TRUE ~ as.character(loc_coded))) %>%
    mutate(ethnic = dplyr::case_when(ethnic == "Mixed/Multiple ethnic groups"~"Other",
                            ethnic == "Asian/Asian British"~"Asian",
                            ethnic == "Black/Black British"~"Black",
                            TRUE ~ as.character(ethnic))) %>%
    mutate(gender = dplyr::case_when(gender == "Non-binary"~"Neither/NA",
                             gender == "Prefer not to say"~"Neither/NA",
                               TRUE ~ as.character(gender)))
}

pipeline_add_factor <- function(data) {
  reg_data <- data %>% mutate(
    age_group = factor(age_group,
      level = c("<30", "30-39", "40-49", "50-59", "60+"),
      labels = c("<30", "30-39", "40-49", "50-59", "60+")),
    gender = factor(gender,
      levels = c("Female", "Male"),
      labels = c("Female", "Male")),
    ethnic = factor(ethnic,
      levels = c("White", "Black", "Asian", "Other"),
      labels = c("White", "Black", "Asian", "Other")),
    c19_diagnosis_category = factor(c19_diagnosis_category),
    job_grouped = factor(job_grouped,
      levels = c("admin",  "allied medical", "domestic", "hca",  "medic",  "nurse",
        "other", "pharmacy", "PT/OT", "Radiographer", NA),
      labels = c("Admin",  "Allied medical", "Domestic", "Healthcare assistants",
        "Medic",  "Nurse", "Other", "Pharmacist", "Occupational and physiotherapists", "Radiographer")),
    loc_coded = factor(loc_coded,
      levels = c("loc_ed", "loc_amu", "loc_critcare", "loc_elderlycare",
        "loc_infectdis", "loc_other", "loc_resp_gsm", "loc_respiratory", NA),
      labels = c("Emergency Department", "Acute Medical Unit", "Critical Care",
        "Geriatric Care",  "Infectious Disease Ward", "Other", "Respiratory Geriatric Ward", "Respiratory Ward")),
    c19zones = factor(c19zones, levels = c("I've not worked in any patient areas",
        "I've worked in green zones most days",
        "I've attended/cared for COVID-19 patients but only after they've been declared no longer infectious",
        "I work in red zones most days but don't attend/care for patients myself",
        "I've attended/cared for patients with COVID-19 in red zones occasionally",
        "I've attended/cared for patients with COVID-19 in red zones most days")))
}

# Conversion from UI to WHO-UI
convert_to_who_ui_s <- function(hero_units) {
  who_units <- vector()
  for (i in seq_len(length(hero_units))) {
    if (hero_units[i] == 0)
      who_units[i] <- 0
    else
      who_units[i] <- 0.0962 + hero_units[i] * 1.0008
  }
  who_units
}

convert_to_who_ui_n <- function(hero_units) {
  who_units <- vector()
  for (i in seq_len(length(hero_units))) {
    if (hero_units[i] == 0)
      who_units[i] <- 0
    else
      who_units[i] <- 0.6501 + hero_units[i] * 0.9348
  }
  who_units
}

pipeline_add_logtitre <- function(data) {
  data_trans <- data %>%
    mutate(IgG_S_titrel2 =
      dplyr::case_when(is.na(IgG_S_titre)~0,
        !is.na(IgG_S_titre)~log2(IgG_S_titre / 100)),
        .after = IgG_S_titre) %>%
    mutate(IgG_N_titrel2 =
      dplyr::case_when(is.na(IgG_N_titre)~0,
        !is.na(IgG_N_titre)~log2(IgG_N_titre / 100)),
          .after = IgG_N_titre) %>%
    mutate(IgG_S_unitl2_hero =
      dplyr::case_when(is.na(IgG_S_unit_hero)~0,
        !is.na(IgG_S_unit_hero)~log2(IgG_S_unit_hero)),
          .after = IgG_S_unit_hero) %>%
    mutate(IgG_N_unitl2_hero =
      dplyr::case_when(is.na(IgG_N_unit_hero)~0,
        !is.na(IgG_N_unit_hero)~log2(IgG_N_unit_hero)),
          .after = IgG_N_unit_hero) %>%
    mutate(IgG_S_unitl2_who =
      convert_to_who_ui_s(IgG_S_unitl2_hero),
        .after = IgG_S_unitl2_hero) %>%
    mutate(IgG_N_unitl2_who =
      convert_to_who_ui_n(IgG_N_unitl2_hero),
        .after = IgG_N_unitl2_hero) %>%
    mutate(IgA_Serum_N_unitl2 =
      dplyr::case_when(is.na(IgA_Serum_N_unit)~0,
        !is.na(IgA_Serum_N_unit)~log2(IgA_Serum_N_unit)),
          .after = IgA_Serum_N_unit) %>%
    mutate(IgA_Serum_S_unitl2 =
      dplyr::case_when(is.na(IgA_Serum_S_unit)~0,
        !is.na(IgA_Serum_S_unit)~log2(IgA_Serum_S_unit)),
          .after = IgA_Serum_S_unit)
}

pipeline_add_date <- function(data) {
  data_trans <- data %>%
    mutate(date_of_sample = lubridate::yday(lubridate::dmy(date_of_sample))) %>%
    mutate(Date_of_PCR = lubridate::yday(lubridate::dmy(Date_of_PCR))) %>%
    mutate(c19_date_of_diagnosis_or_illness =
      lubridate::yday(lubridate::dmy_hms(c19_date_of_diagnosis_or_illness)))
}

pipeline_clean_full <- function(data) {
  data_relabel <- data %>%
    pipeline_relabel_cat %>%
    pipeline_add_factor %>%
    pipeline_add_logtitre %>%
    pipeline_add_date
}

pipeline_new_cols <- function(data) {
  data_trans <- data %>%
    mutate(inf_date =
      dplyr::case_when(!is.na(data$Date_of_PCR)~data$Date_of_PCR,
        !is.na(data$c19_date_of_diagnosis_or_illness)~data$c19_date_of_diagnosis_or_illness),
        .after = c19_date_of_diagnosis_or_illness) %>%
    dplyr::arrange(Sample_ID, date_of_sample) %>%
    mutate(days_since_inf = date_of_sample - inf_date, .after = date_of_sample) %>%
    mutate(symp_pos = dplyr::case_when(is.na(days_since_inf)~0, !is.na(days_since_inf)~1), .after = sero_pos)
}

pipeline_drop_oral <- function(data) {
  data %<>% select(!c(Oralfluid_IgA_spike_of_total, Oralfluid_IgA_NCP_of_total, Oralfluid_IgA_spike_OD,
                      Oralfluid_IgA_spike_units_ngml, Oralfluid_IgA_NCP_OD, Oralfluid_IgA_NCP_units_adjusted_ngml,
                      Oralfluid_IgA_Total_OD, Oralfluid_IgA_Total_units_adjusted_mcgml))
}

pipeline_exclude_drop_outs <- function(data) {
  sample_id_no_follow_up <- data$Sample_ID[is.na(data$sero_pos)]
  data[!data$Sample_ID %in% sample_id_no_follow_up, ]
}

pipeline_exclude_strange_dates <- function(data) {
  data %>% dplyr::filter(is.na(inf_date) | (inf_date > 59 & inf_date < 250))
}

pipeline_get_sero_pos_either <- function(data) {
  ids <- unique(data$HERO_ID)
  bool_sero_pos <- vector()
  data <- data %>% dplyr::arrange(Sample_ID, date_of_sample)
  for (i in seq_len(length(ids))) {
    bool_sero_pos <- c(bool_sero_pos, rep(all(data[data$HERO_ID == ids[i], ][["sero_pos"]] == c(0, 0)), 2))
  }
  data <- data[!bool_sero_pos, ]
  data
}

pipeline_get_sero_revert <- function(data) {
  ids <- unique(data$HERO_ID)
  bool_sero_pos <- vector()
  data <- data %>% dplyr::arrange(Sample_ID, date_of_sample)
  for (i in seq_len(length(ids))) {
    bool_sero_pos <- c(bool_sero_pos, rep(all(data[data$HERO_ID == ids[i], ][["sero_pos"]] == c(1, 0)), 2))
  }
  data <- data[bool_sero_pos, ]
  data
}

pipeline_get_sero_con <- function(data) {
  ids <- unique(data$HERO_ID)
  bool_sero_pos <- vector()
  data <- data %>% dplyr::arrange(Sample_ID, date_of_sample)
  for (i in seq_len(length(ids))) {
    bool_sero_pos <- c(bool_sero_pos, rep(all(data[data$HERO_ID == ids[i], ][["sero_pos"]] == c(0, 1)), 2))
  }
  data <- data[bool_sero_pos, ]
  data
}

pipeline_exclude_sero_convert <- function(data) {
  ids <- unique(data$HERO_ID)
  bool_sero_pos <- vector()
  data <- data %>% dplyr::arrange(Sample_ID, date_of_sample)
  for (i in seq_len(length(ids))) {
    bool_sero_pos <- c(bool_sero_pos, rep(all(data[data$HERO_ID == ids[i], ][["sero_pos"]] == c(1, 0)), 2))
  }
  data <- data[!bool_sero_pos, ]
  data
}

pipeline_get_sero_pos_both <- function(data) {
  ids <- unique(data$HERO_ID)
  bool_sero_pos <- vector()
  data <- data %>% dplyr::arrange(Sample_ID, date_of_sample)
  for (i in seq_len(length(ids))) {
    bool_sero_pos <- c(bool_sero_pos, rep(all(data[data$HERO_ID == ids[i], ][["sero_pos"]] == c(1, 1)), 2))
  }
  data <- data[bool_sero_pos, ]
  data
}

raw_data_dec <- read.csv(system.file("extdata", "HERO_data_clean_21Dec2020.csv", package = "herostudy"))
raw_data_oct <- read.csv(system.file("extdata", "HERO_data_clean_14Oct2020.csv", package = "herostudy"))

raw_data_oct %>% pull(V1_serostatus) %>% table
raw_data_dec %>% pull(V1_serostatus) %>% table

output_rdata <- function() {
  raw_data <- read.csv(system.file("extdata", "HERO_data_clean_21Dec2020.csv", package = "herostudy"))

  custom_pivot <- raw_data %>%
    build_longer_spec(
      cols = contains(c("V1", "V2")),
      names_to = c(".value", "sample_no", "value2"),
      names_pattern = "(.*)V(.)_(.*)"
    )  %>%
    unite(.value, value2, col = ".value", sep = "")
  custom_pivot[31, 2] <- "Serum_IgA_spike_OD"
  custom_pivot[32, 2] <- "Serum_IgA_spike_antibody_units"
  data_clean <- raw_data %>%
    pivot_longer_spec(custom_pivot)  %>%
    dplyr::arrange(Sample_ID) %>%
    pipeline_rename %>%
    pipeline_clean_full %>%
    pipeline_new_cols %>%
    pipeline_drop_oral
  # Exclude samples
  data_clean_sero_base <- data_clean %>%
    pipeline_exclude_drop_outs %>%
    pipeline_exclude_strange_dates
  data_clean_sero <- data_clean_sero_base %>% pipeline_get_sero_pos_both
  data_clean_sero_con <- data_clean_sero_base %>% pipeline_get_sero_con
  data_clean_sero_rev <- data_clean_sero_base %>% pipeline_get_sero_revert
  data_asymp <- data_clean_sero %>% dplyr::filter(is.na(days_since_inf))
  data_symp <- data_clean_sero %>% dplyr::filter(!is.na(days_since_inf))
  data_full <- list(full = data_clean,
                   sero_pos = data_clean_sero,
                   sero_con = data_clean_sero_con,
                   sero_revert = data_clean_sero_rev,
                   symp = data_symp,
                   asymp = data_asymp)
  save(file = here::here("data", "sero_clean.Rdata"), data_full)
  data_full
}

## compare the WHO and HERO universal titres measures and get fitted model
get_linear_rel_ui <- function() {
  y_s <- c(922.74, 568.35, 372.14, 202.75, 116.75, 74.01,
    42.36, 22.15, 12.54, 7.49, 3.64, 1.92)
  y_n <- c(976.32, 532.92, 357.03, 197.64, 136.59, 80.39,
    46.86, 27.35, 14.66, 8.52, 4.90, 3.21)

  x <- vector()
  x[1] <- 1000
  for (i in 1:11) {
     x[i + 1] <- x[i] / 1.75
  }
  m1_s <- lm(y_s ~ x)
  m2_s <- lm(log2(y_s) ~ log2(x))
  summary(m1_s) # 0.9941
  summary(m2_s) # 0.9974, better
  m1_n <- lm(y_n ~ x)
  m2_n <- lm(log2(y_n) ~ log2(x))
  summary(m1_n) # 0.9968
  summary(m2_n) # 0.9984, better
  p1 <- ggplot(data = data.frame(x = x, y = y_S)) +
    geom_point(aes(x = x, y = y), shape = 3, size = 4) +
    geom_abline(intercept = as.numeric(m1_S$coefficients[1]),
      slope = as.numeric(m1_S$coefficients[2]), color = "red") +
    labs(x = "HERO UI, S", y = "WHO UI, S") + theme_bw()
  p2 <- ggplot(data = data.frame(x = x, y = y_N)) +
    geom_point(aes(x = x, y = y), shape = 3, size = 4) +
    geom_abline(intercept = as.numeric(m1_N$coefficients[1]),
      slope = as.numeric(m1_N$coefficients[2]), color = "red") +
    labs(x = "HERO UI, N", y = "WHO UI, N") + theme_bw()
  p3 <- ggplot(data = data.frame(x = log2(x), y = log2(y_S))) +
    geom_point(aes(x = x, y = y), shape = 3, size = 4) +
    geom_abline(intercept = as.numeric(m2_S$coefficients[1]),
      slope = as.numeric(m2_S$coefficients[2]), color = "red") +
    labs(x = "HERO log(UI), S", y = "WHO log(UI), S") + theme_bw()
  p4 <- ggplot(data = data.frame(x = log2(x), y = log2(y_N))) +
    geom_point(aes(x = x, y = y), shape = 3, size = 4) +
    geom_abline(intercept = as.numeric(m2_N$coefficients[1]),
      slope = as.numeric(m2_N$coefficients[2]), color = "red") +
    labs(x = "HERO log(UI), N", y = "WHO log(UI), N") + theme_bw()
  ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
}