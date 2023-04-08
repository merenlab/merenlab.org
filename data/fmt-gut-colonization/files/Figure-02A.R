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

# FUNC PANEL A ------------------------------------------------------------------

plot_detection_vs_prevalence <- function(df_path) {
  
  data_wide <- read_tsv(df_path) %>%
    select(bins, mean_donor_detection, mean_post_detection, CAN_frac_detec)
  
  data <- gather(data_wide, type, mean_detection, mean_donor_detection:mean_post_detection, factor_key=TRUE)
  
  data_donor <- data %>% filter(type == "mean_donor_detection")
  model_donor <- lm(mean_detection ~ CAN_frac_detec, data = data_donor)
  r_squared_donor <- signif(summary(model_donor)$r.squared, 2)
  p_value_donor <- signif(summary(model_donor)$coefficients[,4][2], 2)
  
  data_post <- data %>% filter(type == "mean_post_detection")
  model_post <- lm(mean_detection ~ CAN_frac_detec, data = data_post)
  r_squared_post <- signif(summary(model_post)$r.squared, 2)
  p_value_post <- signif(summary(model_post)$coefficients[,4][2], 2)
  
  plot <- ggplot(data, aes(CAN_frac_detec, mean_detection, color = type, fill = type, linetype = type, shape = type)) +
    geom_smooth(method='lm', formula = y~x, fill="grey60") +
    geom_jitter(height = 0.01, width=0.01, size=0.5) +
    scale_color_manual(values = c("#440154", "#287D8E"), labels = c(paste("Donor: R2=", r_squared_donor, sep=""), paste("Post-FMT: R2=", r_squared_post, sep=""))) +
    scale_fill_manual(values = alpha(c("#440154", "#287D8E"), 1)) +
    scale_shape_manual(values = c(24, 21)) +
    scale_linetype_manual(values=c("solid", "twodash")) +
    xlab("Prevalence in Canada") +
    ylab("Mean Detection") +
    my_theme() +
    guides(fill = FALSE, shape = FALSE, linetype = FALSE)
  
  return(plot)
  
}

# EXECUTION ---------------------------------------------------------------

# Panel A for donor A and donor B
A_DA <- plot_detection_vs_prevalence("mean-detection-vs-prevalence-DA.txt")
A_DB <- plot_detection_vs_prevalence("mean-detection-vs-prevalence-DB.txt")

# make column titles
title1 <- ggdraw() + draw_label("Donor A", fontface='bold')
title2 <- ggdraw() + draw_label("Donor B", fontface='bold')

# make columns
col1 <- plot_grid(title1, A_DA, ncol=1, rel_heights = c(0.1,1))
col2 <- plot_grid(title2, A_DB, ncol=1, rel_heights = c(0.1,1))

# put together
A <- plot_grid(col1, col2, nrow=1)

ggsave("Figure_02A.pdf", A, width=6, height=3, units="in")
