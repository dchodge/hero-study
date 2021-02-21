plot_validation_OD <- function(thres_S, thres_N){
  datafit_OD_val_i <- get_data_val_OD(thres_S, thres_N)
  
  sens_S <- filter(datafit_OD_val_i, aB_protein == "S", Pos_Neg == "Pos") %>% 
    select(sens_ind) %>% summarise(s = mean(sens_ind)*100) %>% round(digits = 2) %>% as.numeric()
  spec_S <- filter(datafit_OD_val_i, aB_protein == "S", Pos_Neg == "Neg") %>% 
    select(sens_ind) %>% summarise(s = mean(sens_ind)*100)  %>% round(digits = 2) %>% as.numeric()
  sens_N <- filter(datafit_OD_val_i, aB_protein == "N", Pos_Neg == "Pos") %>% 
    select(sens_ind) %>% summarise(s = mean(sens_ind)*100)  %>% round(digits = 2) %>% as.numeric()
  spec_N <- filter(datafit_OD_val_i, aB_protein == "N", Pos_Neg == "Neg") %>% 
    select(sens_ind) %>% summarise(s = mean(sens_ind)*100)  %>% round(digits = 2) %>% as.numeric()
  
  ggplot(data = datafit_OD_val_i) + 
    geom_hline(aes(yintercept = OD_thres), linetype = "dashed") + 
    geom_jitter(aes(x = Pos_Neg, y = OD, color = sens_ind), shape = 21, alpha = 0.9, size = 2, width = 0.3) + 
    theme_bw() + labs(x = "Known serostatus of individual", y = expression("OD"[450]), color = "Correctly identified?") + theme(legend.position = "top") +
    facet_grid(cols = vars(aB_protein)) +
    geom_text(aes(x, y , label = paste0("Sens: ", lab)), 
              data=data.frame(x=1, y=4, lab=c(sens_S, sens_N), 
                              aB_protein=c("S","N")), vjust=1) + 
    geom_text(aes(x, y , label = paste0("Spec: ", lab)), 
              data=data.frame(x=2, y=4, lab=c(spec_S, spec_N), 
                              aB_protein=c("S","N")), vjust=1) 
}

plot_validation_OD_HERO <- function(thres_S, thres_N){
  datafit_OD_val_i <- get_data_val_OD(thres_S, thres_N)
  iS1 <- datafit_OD_val_i %>% filter(aB_protein == "S", Pos_Neg == "Pos") %>% pull(sens_ind)
  iN1 <- datafit_OD_val_i %>% filter(aB_protein == "N", Pos_Neg == "Pos") %>% pull(sens_ind)
  iS2 <- datafit_OD_val_i %>% filter(aB_protein == "S", Pos_Neg == "Neg") %>% pull(sens_ind)
  iN2 <- datafit_OD_val_i %>% filter(aB_protein == "N", Pos_Neg == "Neg") %>% pull(sens_ind)
  datafit_OD_val_i$HERO_ind <- c(iS1, iS2|iN2, iN1, iS2|iN2)
  
  sens_S <- filter(datafit_OD_val_i, aB_protein == "S", Pos_Neg == "Pos") %>% 
    select(HERO_ind) %>% summarise(s = mean(HERO_ind)*100) %>% round(digits = 2) %>% as.numeric()
  spec_S <- filter(datafit_OD_val_i, aB_protein == "S", Pos_Neg == "Neg") %>% 
    select(HERO_ind) %>% summarise(s = mean(HERO_ind)*100)  %>% round(digits = 2) %>% as.numeric()
  sens_N <- filter(datafit_OD_val_i, aB_protein == "N", Pos_Neg == "Pos") %>% 
    select(HERO_ind) %>% summarise(s = mean(HERO_ind)*100)  %>% round(digits = 2) %>% as.numeric()
  spec_N <- filter(datafit_OD_val_i, aB_protein == "N", Pos_Neg == "Neg") %>% 
    select(HERO_ind) %>% summarise(s = mean(HERO_ind)*100)  %>% round(digits = 2) %>% as.numeric()
  
  ggplot(data = datafit_OD_val_i) + 
    geom_hline(aes(yintercept = OD_thres), linetype = "dashed") + 
    geom_jitter(aes(x = Pos_Neg, y = OD, color = HERO_ind), shape = 21, alpha = 0.9, size = 2, width = 0.3) + 
    theme_bw() + labs(x = "Known serostatus of individual", y = expression("OD"[450]), color = "Correctly identified?") + theme(legend.position = "top") +
    facet_grid(cols = vars(aB_protein)) +
    geom_text(aes(x, y , label = paste0("Sens: ", lab)), 
              data=data.frame(x=1, y=4, lab=c(sens_S, sens_N), 
                              aB_protein=c("S","N")), vjust=1) + 
    geom_text(aes(x, y , label = paste0("Spec: ", lab)), 
              data=data.frame(x=2, y=4, lab=c(spec_S, spec_N), 
                              aB_protein=c("S","N")), vjust=1) 
}

plot_OD_sens_curve <- function(){
  
  df_label <- data.frame(x = c(50, 50, 50, 50), y = c(0.3, 2.1, 0.3, 3), protein = c("S","S","N","N"))
  
  ggplot(data = datafit_OD_curve, aes(y = OD)) +
    geom_point(aes(x = Value, color = Type) , shape = 4) +
    labs(y = expression("OD"[450]), x = "%", fill = "") + 
    scale_fill_identity(guide = 'legend', labels = c('Specificity', 'Sensivitiy')) + 
    theme_bw() + theme(legend.position="bottom") + facet_grid(cols = vars(protein)) + 
    geom_text(data = df_label, aes(x = x, y = y), label = c(expression("f"["sp"]^"S"),
                                                            expression("f"["se"]^"S"),
                                                            expression("f"["sp"]^"N"),
                                                            expression("f"["se"]^"N")))
}



plot_ROC <- function(){
  datafit_ROC <- datafit_OD_curve %>% select(OD, Value, Type, protein) %>% 
    tidyr::spread(Type, Value) %>% mutate(FNR = (100-Sensitivity)/100, FPR = (100-Specificity)/100)
  
  lol <- data.frame()
  for (i in seq(0, 2, 0.2))
  {
    lol <- rbind(lol, filter(datafit_ROC, protein=="N")[which.min(abs(filter(datafit_ROC, protein=="N")$OD-i)),])
    lol <- rbind(lol, filter(datafit_ROC, protein=="S")[which.min(abs(filter(datafit_ROC, protein=="S")$OD-i)),])
  }
  lol$OD_label <-  seq(0, 2, 0.2) %>% purrr::map(~rep(.x, 2)) %>% unlist
  
  ggplot(data = datafit_ROC, aes(x = FPR, y = Sensitivity/100)) +
    geom_line(shape = 4) + geom_abline() + 
    theme_bw() + theme(legend.position="bottom") + 
    geom_point(data = lol, aes(y = Sensitivity/100, x = FPR), color = palette(rainbow(22))) +
    geom_text(data = lol, aes(y = Sensitivity/100 - 0.024, x = FPR + 0.04, label = OD_label), color = 'red', size = 2.5) +
    facet_grid(cols = vars(protein)) 
}


get_dataplt_sens <- function(sensvals){
  
  prev_df <- vector()
  sens_df <- vector()
  df_best <- data.frame()
  
  for (p in c("S", "N")) {
    hero_OD <- filter(datafit_hero, sero_pos == 1) %>% na.omit %>% select(age_group, paste0("IgG_", p, "_OD"))
    datafit_OD_Sens <- datafit_OD_curve %>% filter(Type == "Sensitivity", protein == p)
    for (s in sensvals){
      probs1 <- rbeta(100, 190*s/100, 190*(100 - s)/100)
      df_OD <- data.frame()
      for (prob in probs1){
        OD_value <- datafit_OD_Sens[which.min(abs(datafit_OD_Sens$Value - prob*100)),]
        df_OD <- rbind(df_OD, OD_value)
      }
      OD_sens <- df_OD$OD
      for (i in OD_sens){
        prev_i <- hero_OD %>% mutate(ind = (eval(parse(text = paste0("IgG_", p, "_OD"))) > i)) %>% 
          group_by(age_group) %>% 
          summarise(pos = sum(ind), n = n()) %>% mutate(sens = pos/n)  %>%
          pull(sens)
        prev_df <- c(prev_df, prev_i)
      }
      sens_df <- c(sens_df, rep(s, 5*100))
    }
    
    df_best_pro <- data.frame(sens = sens_df, implied_sens = prev_df*100, 
                              age_group = c("<30", "30-39", "40-49", "50-59", "60+"),
                              protein = p)
    
    
    df_best <- rbind(df_best_pro, df_best)
  }
  
  dataplt_sens <- df_best %>% group_by(age_group, sens, protein) %>% summarise(implied_sens = mean(implied_sens)) 
  prop_age <- (table(hero_OD$age_group)/length(hero_OD$age_group)) %>% as.numeric
  dataplt_sens_all <- dataplt_sens %>% group_by(sens, protein) %>% summarise(mean = weighted.mean(implied_sens, prop_age))
  dataplt_sens <- list(age = dataplt_sens, all = dataplt_sens_all)
  
}


plot_het_sens <- function(){
  
  ggplot(data = dataplt_sens$age) + 
    geom_abline(linetype = "dashed") + 
    geom_line(data = dataplt_sens$all, aes(x = sens, y = mean), color = "black", size = 1, alpha = 0.6)   + 
    geom_line(aes(x = sens, y = prev, color = age_group), size = 1, alpha = 0.9) + theme_bw()  + 
    labs(x = "Sensitivity of controls", y = "Sensitivity of HERO samples") + 
    facet_grid(rows = vars(protein))
}

get_dataplt_spec <- function(sensvals){
  
  prev_df <- vector()
  sens_df <- vector()
  df_best <- data.frame()
  for (p in c("S", "N")){
    hero_OD <- filter(datafit_hero, sero_pos == 0) %>% na.omit %>% 
      select(age_group, paste0("IgG_", p, "_OD"))
    datafit_OD_Sens <- datafit_OD_curve %>% filter(Type == "Specificity", protein == p)
    for (s in sensvals){
      probs1 <- rbeta(100, 690*s/100, 690*(100 - s)/100)
      df_OD <- data.frame()
      for (prob in probs1){
        OD_value <- datafit_OD_Sens[which.min(abs(datafit_OD_Sens$Value - prob*100)),]
        df_OD <- rbind(df_OD, OD_value)
      }
      OD_sens <- df_OD$OD
      for (i in OD_sens){
        prev_i <- hero_OD %>% mutate(ind = (eval(parse(text = paste0("IgG_", p, "_OD"))) < i)) %>% group_by(age_group) %>% 
          summarise(pos = sum(ind), n = n()) %>% mutate(sens = pos/n)  %>%
          pull(sens)
        prev_df <- c(prev_df, prev_i)
      }
      sens_df <- c(sens_df, rep(s, 5*100))
    }
    
    df_best_pro <- data.frame(sens = sens_df, prev = prev_df*100, 
                              age_group = c("<30", "30-39", "40-49", "50-59", "60+"),
                              protein = p)
    df_best <- rbind(df_best_pro, df_best)
  }
  dataplt_spec <- df_best %>% group_by(age_group, sens, protein) %>% summarise(prev = mean(prev)) 
  prop_age <- (table(hero_OD$age_group)/length(hero_OD$age_group)) %>% as.numeric
  dataplt_spec_all <- dataplt_spec %>% group_by(sens, protein) %>% summarise(mean = weighted.mean(prev, prop_age))
  dataplt_spec <- list(age = dataplt_spec, all = dataplt_spec_all)
}

plot_het_spec <- function(){
  ggplot(data = dataplt_spec$age) + 
    geom_abline(linetype = "dashed") + 
    geom_line(data = dataplt_spec$all, aes(x = sens, y = mean), color = "black", size = 1, alpha = 0.6)   + 
    geom_line(aes(x = sens, y = prev, color = age_group), size = 1, alpha = 0.9) + theme_bw()  + 
    labs(x = "Specificity of controls", y = "Specificity of HERO samples") + 
    facet_grid(rows = vars(protein))
}

pltfigure3 <- function(height, width){
  p1 <- plot_het_sens()
  p2 <- plot_het_spec()
  ggpubr::ggarrange(p1, p2, ncol = 2, common.legend = TRUE)
  ggsave(filename = here::here("outputs","figs", paste0("fig3", ".pdf")), height = height, width = width)
}


plot_ROC_age <- function(){
  datafit_ROC <- datafit_OD_curve %>% select(OD, Value, Type, protein) %>% 
    tidyr::spread(Type, Value) %>% mutate(FNR = (100-Sensitivity)/100, FPR = (100-Specificity)/100)
  
  datafit_raw <- datafit_hero %>% select(IgG_S_OD, IgG_N_OD, age_group, sero_pos)
  base_sens <- datafit_raw %>% filter(sero_pos == 1) %>% group_by(age_group) %>% 
    summarise(freq = sum(sero_pos)) %>% pull(freq)
  base_spec <- datafit_raw %>% filter(sero_pos == 0) %>% group_by(age_group) %>% 
    mutate(sero_pos = 1) %>% summarise(freq = sum(sero_pos)) %>% pull(freq)
  
  datafit_full <- data.frame()
  for (i in seq(0, 3, 0.01)){
    ODsum_sens <-  datafit_raw %>% filter(sero_pos == 1) %>% mutate(sero_pos_pred = (i < IgG_S_OD)) %>%
      group_by(age_group) %>% summarise(freq = sum(sero_pos_pred)) %>% mutate(sens = freq/base_sens, OD = i, protein = "S") %>% 
      select(age_group, sens, protein, OD)
    ODsum_spec <- datafit_raw %>% filter(sero_pos == 0) %>% mutate(sero_pos_pred = (i > IgG_S_OD)) %>%
      group_by(age_group) %>% summarise(freq = sum(sero_pos_pred)) %>% mutate(spec = freq/base_spec, OD = i, protein = "S") %>% 
      select(age_group, spec, protein, OD)
    ODsumS <- merge(ODsum_sens, ODsum_spec, by = c("age_group", "OD", "protein"))
    
    ODsum_sens <-  datafit_raw %>% filter(sero_pos == 1) %>% mutate(sero_pos_pred = (i < IgG_N_OD)) %>%
      group_by(age_group) %>% summarise(freq = sum(sero_pos_pred)) %>% mutate(sens = freq/base_sens, OD = i, protein = "N") %>% 
      select(age_group, sens, protein, OD)
    ODsum_spec <- datafit_raw %>% filter(sero_pos == 0) %>% mutate(sero_pos_pred = (i > IgG_N_OD)) %>%
      group_by(age_group) %>% summarise(freq = sum(sero_pos_pred)) %>% mutate(spec = freq/base_spec, OD = i, protein = "N") %>% 
      select(age_group, spec, protein, OD)
    ODsumN <- merge(ODsum_sens, ODsum_spec, by = c("age_group", "OD", "protein"))
    
    datafit_full <- rbind(datafit_full, ODsumS, ODsumN)
  }
  datafit_full %<>% select(age_group, sens, spec, protein, OD) %>% unique %>% mutate(FPR = 1-spec) %>% arrange(sens)
  
  lol <- data.frame()
  for (i in seq(0, 2, 0.2))
  {
    lol1 <- datafit_full %>% filter(OD == i)
    lol <- rbind(lol, lol1)
  }
  
  ggplot(data = datafit_full, aes(x = FPR, y = sens)) +
    geom_line(size = 0.5, alpha = 0.8) + 
    theme_bw() + theme(legend.position="bottom") + 
    labs(x = "FPR", y = "Sensitivity") + facet_grid(rows = vars(protein), cols = vars(age_group)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    geom_point(data = lol, aes(y = sens, x = FPR), color = (palette(rainbow(9)) %>% purrr::map(~rep(.x, 10)) %>% unlist)) +
    geom_text(data = lol, aes(y = sens - 0.02, x = FPR + 0.02, label = OD), color = (palette(rainbow(9)) %>% purrr::map(~rep(.x, 10)) %>% unlist), size = 2.5, position=position_jitter(width=0.003,height=0.001)) +
    xlim(0, 0.2) + ylim(0.25, 1)
  
}