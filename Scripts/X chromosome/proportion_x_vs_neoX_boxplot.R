library(tidyverse)

# Inputs
x_dir <- "Genes_on_X_noNeoX/bed"
neo_dir <- "Genes_on_neoX/bed"
stats_file <- "genome_stats_by_species.csv"

count_genes <- function(bed_dir, label) {

  bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

  if (length(bed_files) == 0) {
    stop(paste("No .bed files found in:", bed_dir))
  }

  tibble(
    genome = basename(bed_files),
    x_genes = sapply(
      bed_files,
      function(f) {
        nrow(read.table(
          f,
          sep = "\t",
          header = FALSE,
          stringsAsFactors = FALSE,
          quote = "",
          comment.char = ""
        ))
      }
    ),
    group = label
  ) %>%
    mutate(genome = sub("\\.bed$", "", genome))
}

#Count genes
x_counts <- count_genes(x_dir, "Ancestral X")
neo_counts <- count_genes(neo_dir, "neoX")

gene_counts <- bind_rows(x_counts, neo_counts)

#Read genome stats
stats <- read_csv(stats_file, show_col_types = FALSE) %>%
  transmute(
    genome = sub("\\.gff3$", "", species_file),
    genes
  )

#Join + calculate proportion
plot_df <- gene_counts %>%
  left_join(stats, by = "genome") %>%
  mutate(proportion = 100 * x_genes / genes)

#Check mismatches
missing <- plot_df %>% filter(is.na(genes))
if (nrow(missing) > 0) {
  cat("These genomes still did not match:\n")
  print(missing$genome)
}

plot_df$group <- factor(plot_df$group, levels = c("Ancestral X", "neoX"))
plot_df$xpos <- ifelse(plot_df$group == "Ancestral X", 1, 1.6)

# Plot
p <- ggplot(plot_df, aes(x = xpos, y = proportion, group = group)) +

  geom_boxplot(width = 0.25, outlier.shape = NA) +

  geom_jitter(width = 0.05, size = 2, alpha = 0.7) +

  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.245,
    color = "red",
    linetype = "solid",
    linewidth = 0.15
  ) +

  scale_x_continuous(
    breaks = c(1, 1.6),
    labels = c("Ancestral X", "neoX")
  ) +

  labs(
    x = "",
    y = "Proportion of genes (%)",
    title = "Proportion of X-linked genes"
  ) +

  theme_classic(base_size = 14)

ggsave("X_vs_neoX_gene_proportion_boxplot.pdf", p, width = 6, height = 6)
ggsave("X_vs_neoX_gene_proportion_boxplot.png", p, width = 6, height = 6, dpi = 300)
