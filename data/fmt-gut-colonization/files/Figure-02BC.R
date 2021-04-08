<<<<<<< HEAD
#!/usr/bin/env Rscript
=======
>>>>>>> 05b25422ae07009bc26f00953a8f33577d51e51c
library(tidyverse)
library(pROC)
library(cowplot)
library(ggsignif)


setwd("/Users/andrea/github/web/data/fmt-gut-colonization/files")

# THEME -------------------------------------------------------------------

my_theme <- function () {
  theme_bw() + (theme(plot.caption = element_blank(), 
                      plot.title = element_blank(),
                      panel.border = element_rect(colour = "black", fill=NA, size=1),
                      text = element_text(colour = "black", family="Helvetica"),
                      axis.text = element_text(color="black", size = 7),
                      legend.key.size = unit(1, 'lines'),
                      legend.text = element_text(size=9),
                      legend.title = element_blank(),
                      legend.margin=margin(c(0,0,4,4)),
                      legend.background = element_rect(fill=alpha("white", 0.5)),
                      legend.justification = c(1,0),
                      legend.position = c(1, 0),
                      axis.title = element_text(size=10, face="bold"),
                      plot.margin = unit(c(8,8,8,8), "pt")))
}

# FUNC PANEL B & PANEL C ------------------------------------------------------------------

plot_logit_and_auc <- function (df_path) {
  
  # Reading df & releveling outcome variable:
  df <- read_tsv(df_path) 
  df <- df %>% 
    mutate(outcome = fct_relevel(outcome, "no_colonization", "colonization")) %>%
    mutate(transplant_mean_cov_Q2Q3 = log(transplant_mean_cov_Q2Q3, base = 10))
  
  # LOGISTIC REGRESSION MODELS
  
  # mct = mean coverage in transplant sample
  # cfd = fraction detected in Canada 
  logitMod_mct <- glm(outcome ~ transplant_mean_cov_Q2Q3, data = df, family = binomial("logit"))
  logitMod_cfd <- glm(outcome ~ CAN_frac_detected, data = df, family = binomial("logit"))
  logitMod_all <- glm(outcome ~ transplant_mean_cov_Q2Q3 + CAN_frac_detected, data = df, family = binomial("logit"))
  
  p_value_mct <- signif(summary(logitMod_mct)$coefficients[2,4], 2)
  p_value_cfd <- signif(summary(logitMod_cfd)$coefficients[2,4], 2)
  
  # BOXPLOTS
  
  # Make informative labels with n values:
  colo_label <- paste("COLO\nn =", sum(df$outcome == "colonization"), sep = " ")
  no_colo_label <- paste("NO COLO\nn =", sum(df$outcome == "no_colonization"), sep = " ")
  
  if (p_value_mct >= 0.05) {
    sig <- "n.s."
  } else if (p_value_mct >= 0.01) {
    sig <- "*"
  } else if (p_value_mct >= 0.001) {
    sig <- "**"
  } else if (p_value_mct < 0.001) {
    sig <- "***"
  }
  
  # Box plot for mean coverage Q2Q3 in transplant sample:
  bp_mct <- ggplot(df, aes(x = transplant_mean_cov_Q2Q3, y = outcome)) +
    geom_boxplot(width = 0.75, size = 0.5, color = alpha("black", 1), alpha = 0, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.01, height = 0.3), size = 0.5, fill = alpha("black", 1), color = alpha("black", 1), shape = 21) +
    labs(title = "Dose", caption = paste("p-value:", p_value_mct)) +
    xlab("log10(Coverage in Transplant)") +
    scale_y_discrete(labels = c(no_colo_label, colo_label)) +
    my_theme() +
    theme(axis.title.y = element_blank()) +
    geom_segment(aes(x=min(df$transplant_mean_cov_Q2Q3)-0.05, y=1, xend=min(df$transplant_mean_cov_Q2Q3)-0.05, yend=2), size=.5) +
    annotate(geom="text", x=min(df$transplant_mean_cov_Q2Q3)-0.05, y=1.5, size=3.5, angle=90, label=sig)
  
  if (p_value_cfd >= 0.05) {
    sig <- "n.s."
  } else if (p_value_cfd >= 0.01) {
    sig <- "*"
  } else if (p_value_cfd >= 0.001) {
    sig <- "**"
  } else if (p_value_cfd < 0.001) {
    sig <- "***"
  }
  
  # Box plot for fraction detected in Canada:
  bp_cfd <- ggplot(df, aes(x = CAN_frac_detected, y = outcome)) +
    geom_boxplot(width = 0.75, size = 0.5, color = alpha("black", 1), alpha = 0, outlier.shape = NA) +
    geom_point(position = position_jitter(width = 0.01, height = 0.3), size = 0.5, fill = alpha("black", 1), color = alpha("black", 1), shape = 21) +
    labs(title = "Prevalence", caption = paste("p-value:", p_value_cfd)) +
    xlab("Prevalence in Canada") +
    scale_y_discrete(labels = c(no_colo_label, colo_label)) +
    my_theme() +
    theme(axis.title.y = element_blank()) +
    geom_segment(aes(x=min(df$CAN_frac_detected)-0.05, y=1, xend=min(df$CAN_frac_detected)-0.05, yend=2), size=.5) +
    annotate(geom="text", x=min(df$CAN_frac_detected)-0.05, y=1.5, size=3.5, angle=90, label=sig)
  
  # ROC CURVES
  
  predictions_mct <- predict(logitMod_mct, df["transplant_mean_cov_Q2Q3"], outcome="response")
  roc_obj_mct <- roc(df$outcome, predictions_mct)
  auc_mct <- signif(roc_obj_mct$auc, digits = 2)
  roc_df_mct <- data.frame(sensitivity = roc_obj_mct$sensitivities, 
                           specificity = roc_obj_mct$specificities, 
                           threshold = roc_obj_mct$thresholds,
                           vars = paste("Dose\nAUC=", auc_mct, sep=""))
  
  predictions_cfd <- predict(logitMod_cfd, df["CAN_frac_detected"], outcome="response")
  roc_obj_cfd <- roc(df$outcome, predictions_cfd)
  auc_cfd <- signif(roc_obj_cfd$auc, digits = 2)
  roc_df_cfd <- data.frame(sensitivity = roc_obj_cfd$sensitivities, 
                           specificity = roc_obj_cfd$specificities, 
                           threshold = roc_obj_cfd$thresholds,
                           vars = paste("Prevalence\nAUC=", auc_cfd, sep=""))
  
  predictions_all <- predict(logitMod_all, df[c("transplant_mean_cov_Q2Q3", "CAN_frac_detected")], outcome="response")
  roc_obj_all <- roc(df$outcome, predictions_all)
  auc_all <- signif(roc_obj_all$auc, digits = 2)
  roc_df_all <- data.frame(sensitivity = roc_obj_all$sensitivities, 
                           specificity = roc_obj_all$specificities, 
                           threshold = roc_obj_all$thresholds,
                           vars = paste("Dose + Prevalence\nAUC=", auc_all, sep=""))
  
  roc_df <- bind_rows(roc_df_mct, roc_df_cfd, roc_df_all)
  
  p_roc <- ggplot(roc_df, aes(x = 1 - specificity, y = sensitivity, group = vars)) +
    geom_abline(color="grey90") +
    geom_path(size = 1, aes(color = vars, linetype = vars)) +
    scale_color_manual(values = c("#FFC161", "#C54C4C", "#386E97")) +
    scale_linetype_manual(values=c("solid", "dotdash", "dotted")) +
    xlab("1 - Specificity") +
    ylab("Sensitivity") +
    labs(title = "ROC Curves of Colonization Predictors", linetype = "Predictor Variable(s)", color = "Predictor Variable(s)") +
    my_theme() +
    theme(legend.background = element_blank())
  
  return(list(bp_mct, bp_cfd, p_roc))
  
}

# EXECUTION ---------------------------------------------------------------

# Panel B and panel C for donor A and donor B
BC_DA <- plot_logit_and_auc(df_path = "summary-DA.txt")
BC_DB <- plot_logit_and_auc(df_path = "summary-DB.txt")

B1_DA <- BC_DA[[1]]
B2_DA <- BC_DA[[2]]
C_DA <- BC_DA[[3]]

B1_DB <- BC_DB[[1]]
B2_DB <- BC_DB[[2]]
C_DB <- BC_DB[[3]]


# Panel A for donor A and donor B
# Panel B and panel C for donor A and donor B
align_first_col <- align_plots(B1_DA, B2_DA, C_DA, align='v', axis='l')
align_second_col <- align_plots(B1_DB, B2_DB, C_DB, align='v', axis='l')

# squish boxplots
B_DA <- plot_grid(align_first_col[[1]], align_first_col[[2]], ncol=1)
B_DB <- plot_grid(align_second_col[[1]], align_second_col[[2]], ncol=1)

# make column titles
title1 <- ggdraw() + draw_label("Donor A", fontface='bold')
title2 <- ggdraw() + draw_label("Donor B", fontface='bold')

# make columns
col1 <- plot_grid(title1, B_DA, align_first_col[[3]], ncol=1, rel_heights = c(0.1,1,1))
col2 <- plot_grid(title2, B_DB, align_second_col[[3]], ncol=1, rel_heights = c(0.1,1,1))
# put together
BC <- plot_grid(col1, col2, nrow=1)

<<<<<<< HEAD
ggsave("Figure-02BC.pdf", BC, width=6, height=6, units="in", dpi=300)
=======
ggsave("Figure-02BC.pdf", BC, width=6, height=6, units="in", dpi=300)
>>>>>>> 05b25422ae07009bc26f00953a8f33577d51e51c
