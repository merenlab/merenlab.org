library(tidyverse)
library(cowplot) 

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


# FUNC PANEL D ------------------------------------------------------------

plot_dose_vs_prevalence <- function (df_path) {
  
  # Reading df, taking log of coverage, releveling outcome variable:
  df <- read_tsv(df_path)
  df <- df %>% 
    mutate(transplant_mean_cov_Q2Q3 = log(transplant_mean_cov_Q2Q3, base = 10)) %>%
    mutate(outcome = fct_relevel(outcome, "no_colonization", "colonization"))
  
  # Making linear model and getting statistics:
  model <- lm(transplant_mean_cov_Q2Q3 ~ CAN_frac_detected, data = df)
  r_squared <- signif(summary(model)$r.squared, 2)
  p_value <- signif(summary(model)$coefficients[,4][2], 2)
  
  # Making plots:  
  p <- ggplot(df, aes(CAN_frac_detected, transplant_mean_cov_Q2Q3)) +
    geom_smooth(method = 'lm', formula = y~x, color = "black") +
    geom_jitter(shape = 21, height = 0.01, width=0.01, fill = alpha("black", 1), color = alpha("black", 1), size = 0.5) +
    xlab("Prevalence in Canada") +
    ylab("log10(Coverage in Transplant)") +
    my_theme() +
    annotate(geom="text", x=1, y=-1, hjust=1, vjust=0, size=3.5, label=paste("R2=", r_squared, " p=", p_value, sep=""))
  
  return(p)
  
}


# EXECUTION ---------------------------------------------------------------

# Panel D for donor A and donor B
D_DA <- plot_dose_vs_prevalence(df_path = "summary-DA.txt")
D_DB <- plot_dose_vs_prevalence(df_path = "summary-DB.txt")

D <- plot_grid(D_DA, D_DB, nrow=1)

ggsave("Figure_02D.pdf", D, width=6, height=3, units="in")