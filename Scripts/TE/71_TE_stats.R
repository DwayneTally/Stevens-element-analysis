library(tidyverse)
library(broom)

df <- read_csv("71_repeat_summary.csv")

df_long <- df %>%
  pivot_longer(
    cols = c("DNA", "Rolling Circle", "Penelope", "LINE", "SINE", "LTR", "Simple Repeat (Other)", "Unclassified"),
    names_to = "TE_class",
    values_to = "percent"
  )

model_stats <- df_long %>%
  group_by(TE_class) %>%
  group_modify(~{
    model <- lm(percent ~ log10(`Genome Size`), data = .x)

    tidy_mod <- tidy(model)
    glance_mod <- glance(model)

    slope_row <- tidy_mod %>%
      filter(term != "(Intercept)") %>%
      slice(1)

    tibble(
      slope = slope_row$estimate,
      p_value = slope_row$p.value,
      r_squared = glance_mod$r.squared
    )
  }) %>%
  ungroup() %>%
  mutate(
    label = paste0(
      "slope = ", signif(slope, 3), "\n",
      "R² = ", signif(r_squared, 3), "\n",
      "p = ", signif(p_value, 3)
    )
  )

p <- ggplot(df_long, aes(x = `Genome Size`, y = percent)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ TE_class) +
  geom_text(
    data = model_stats,
    aes(x = Inf, y = Inf, label = label),
    hjust = 1.1,
    vjust = 1.1,
    inherit.aes = FALSE,
    size = 3.5
  ) +
  scale_x_log10() +
  labs(
    x = "Genome Size (bp, log10 scale)",
    y = "Percentage of TE Class",
    title = "TE Composition vs Genome Size"
  ) +
  theme_classic(base_size = 14)

print(p)

ggsave("TE_vs_genome_size_with_stats.pdf", p, width = 12, height = 8, dpi = 600)

# optional: print the stats table too
print(model_stats)
write_csv(model_stats, "TE_vs_genome_size_model_stats.csv")
