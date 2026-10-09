library(tidyverse)

bacteria <- read_csv("bacteria.csv", show_col_types = FALSE) %>%
  mutate(present = if_else(y == "y", 1, 0),
         trt = factor(trt, levels = c("placebo", "drug", "drug+")))

# proportion of positive tests per treatment
summary_table <- bacteria %>%
  group_by(trt) %>%
  summarise(n_tests = n(),
            percent_positive = round(100 * mean(present), 1))
print(summary_table)

# logistic regression: does treatment reduce H. influenzae presence?
model <- glm(present ~ trt + week, family = binomial, data = bacteria)
print(summary(model)$coefficients)

p <- ggplot(summary_table, aes(x = trt, y = percent_positive, fill = trt)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = paste0(percent_positive, "%")), vjust = -0.5) +
  scale_fill_manual(values = c("grey60", "#1b9e77", "#7570b3")) +
  labs(x = NULL, y = "Tests positive for H. influenzae (%)",
       title = "Drug treatment reduces H. influenzae carriage") +
  ylim(0, 100) +
  theme_minimal() + theme(legend.position = "none")
ggsave("barplot.png", p, width = 5, height = 4, dpi = 150)
