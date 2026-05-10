#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(dplyr)
})

df <- read_csv("genome_stats_by_species.csv", show_col_types = FALSE)

# Normalize Tenebrionidea naming
df$group <- recode(df$group,
  "Tenebrionidea" = "Tenebrionoidea",
  "Tenebrionidae" = "Tenebrionoidea"
)

# Explicit color palette (superfamily -> color) IN DESIRED ORDER
sf_colors <- c(
  "Caraboidea"      = "#d5a6bd",
  "Buprestoidea"    = "#783f04",
  "Elateroidea"     = "#6aa84f",
  "Staphylinoidea"  = "#aaaaaa",
  "Scarabaeoidea"   = "#f6b26b",
  "Cleroidea"       = "#cfe2f3",
  "Coccinelloidea"  = "#f1c232",
  "Tenebrionoidea"  = "#3d85c6",
  "Cucujoidea"      = "#8e7cc3",
  "Curculionoidea"  = "#e69138",
  "Chrysomeloidea"  = "#b6d7a8"
)

# Keep only groups with >= 2 species
df2 <- df %>%
  group_by(group) %>%
  filter(n() >= 2) %>%
  ungroup()

# Drop any groups not in the color map (safety)
df2 <- df2 %>% filter(group %in% names(sf_colors))

# Force x-axis order to match sf_colors order
df2$group <- factor(df2$group, levels = names(sf_colors))

# Long format for plotting (factor levels carry through)
long <- df2 %>%
  select(group, species_file, chromosomes, genome_size_mb, gc_percent, genes) %>%
  pivot_longer(
    cols = c(chromosomes, genome_size_mb, gc_percent, genes),
    names_to = "metric",
    values_to = "value"
  )

metric_labs <- c(
  chromosomes    = "Chromosomes (N)",
  genome_size_mb = "Genome size (Mbp)",
  gc_percent     = "GC (%)",
  genes          = "Genes "
)

p <- ggplot(long, aes(x = group, y = value)) +
  geom_boxplot(
    aes(fill = group),
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(color = group),
    width = 0.18,
    height = 0,
    size = 1.3,
    alpha = 0.8
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    labeller = as_labeller(metric_labs)
  ) +
  scale_fill_manual(values = sf_colors, drop = FALSE) +
  scale_color_manual(values = sf_colors, drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  )

ggsave(
  "genome_stats_boxplots_by_superfamily_colored.pdf",
  p,
  width = 14,
  height = 8
)
ggsave(
  "genome_stats_boxplots_by_superfamily_colored.svg",
  p,
  width = 14,
  height = 8,
  dpi = 300
)

message("Wrote genome_stats_boxplots_by_superfamily_colored.pdf and .svg")

