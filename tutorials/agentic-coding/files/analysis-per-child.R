library(tidyverse)

bacteria <- read_csv("bacteria.csv", show_col_types = FALSE) %>%
  mutate(present = if_else(y == "y", 1, 0),
         trt = factor(trt, levels = c("placebo", "drug", "drug+")))

# how many rows, and how many children?
cat("rows:", nrow(bacteria), "  children:", n_distinct(bacteria$ID), "\n")

# how many tests per child?
print(table(table(bacteria$ID)))

# one number per child: the fraction of their tests that were positive
per_child <- bacteria %>%
  group_by(ID, trt) %>%
  summarise(fraction_positive = mean(present), n_tests = n(), .groups = "drop")
print(count(per_child, trt))

# drug (either arm) vs placebo, one value per child
per_child <- per_child %>% mutate(any_drug = trt != "placebo")
print(wilcox.test(fraction_positive ~ any_drug, data = per_child, exact = FALSE))

p <- ggplot(bacteria, aes(x = factor(week), y = fct_reorder(ID, present, .fun = mean), fill = y)) +
  geom_tile(color = "white") +
  facet_grid(trt ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = c(n = "grey85", y = "#d95f02"), labels = c("absent", "present")) +
  labs(x = "Week", y = "Child", fill = "H. influenzae") +
  theme_minimal() + theme(axis.text.y = element_blank(), panel.grid = element_blank())
ggsave("per-child-tiles.png", p, width = 5, height = 7, dpi = 150)

set.seed(42)
p2 <- ggplot(per_child, aes(x = trt, y = fraction_positive)) +
  geom_boxplot(outlier.shape = NA, width = 0.5, color = "grey50") +
  geom_jitter(aes(size = n_tests), width = 0.15, height = 0, alpha = 0.6) +
  scale_size_continuous(range = c(1.5, 4), breaks = 2:5) +
  labs(x = NULL, y = "Fraction of a child's tests that were positive",
       size = "Tests per child") +
  theme_minimal()
ggsave("per-child-dots.png", p2, width = 5, height = 4, dpi = 150)
