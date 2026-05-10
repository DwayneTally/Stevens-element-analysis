library(tidyverse)
library(ggrepel)

bed_dir <- "/N/project/Bracewell_fly/Dwayne/helixer_results/Genes_on_X_noNeosex/bed"
map_file <- "/N/project/Bracewell_fly/Dwayne/helixer_results/genome_superfamily.tsv"

bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

if (length(bed_files) == 0) {
  stop("No .bed files found in: ", bed_dir)
}

# Count genes per BED file
gene_counts <- tibble(
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
  )
) %>%
  mutate(genome = sub("\\.bed$", "", genome))

# Read genome -> superfamily map
map_df <- read_tsv(map_file, show_col_types = FALSE)

# Merge
plot_df <- gene_counts %>%
  left_join(map_df, by = "genome")

# Check for missing superfamily assignments
missing_sf <- plot_df %>% filter(is.na(superfamily))
if (nrow(missing_sf) > 0) {
  cat("These genomes are missing a superfamily assignment:\n")
  print(missing_sf$genome)
}

# Define outliers within each superfamily using 1.5*IQR
plot_df <- plot_df %>%
  group_by(superfamily) %>%
  mutate(
    Q1 = quantile(x_genes, 0.25, na.rm = TRUE),
    Q3 = quantile(x_genes, 0.75, na.rm = TRUE),
    IQRv = IQR(x_genes, na.rm = TRUE),
    lower = Q1 - 1.5 * IQRv,
    upper = Q3 + 1.5 * IQRv,
    is_outlier = x_genes < lower | x_genes > upper
  ) %>%
  ungroup()

# Optional: order superfamilies by median x_genes
sf_order <- plot_df %>%
  group_by(superfamily) %>%
  summarize(med = median(x_genes, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(med)) %>%
  pull(superfamily)

plot_df$superfamily <- factor(plot_df$superfamily, levels = sf_order)

# Plot
p <- ggplot(plot_df, aes(x = superfamily, y = x_genes, fill = superfamily)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.8, width = 0.7) +
  geom_jitter(aes(color = superfamily), width = 0.15, size = 2, alpha = 0.75, show.legend = FALSE) +
  geom_text_repel(
    data = plot_df %>% filter(is_outlier),
    aes(label = genome),
    size = 3,
    max.overlaps = 50,
    box.padding = 0.35,
    point.padding = 0.2,
    segment.color = "grey40",
    show.legend = FALSE
  ) +
  labs(
    x = "Superfamily",
    y = "Number of genes",
    title = "X-linked gene counts"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  )

print(p)

ggsave("X_gene_counts_by_superfamily_boxplot.pdf", p, width = 10, height = 7)
ggsave("X_gene_counts_by_superfamily_boxplot.png", p, width = 10, height = 7, dpi = 300)

# Also save the table with outlier calls
write_tsv(plot_df, "X_gene_counts_with_superfamily_and_outliers.tsv")
