library(tidyverse)

# Directories
x_dir <- "Genes_on_X_noNeoX/bed"
neo_dir <- "Genes_on_neoX/bed"

count_genes <- function(bed_dir, label) {

  bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

  if (length(bed_files) == 0) {
    stop(paste("No .bed files found in:", bed_dir))
  }

  df <- tibble(
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
  )

  df$genome <- sub("\\.bed$", "", df$genome)
  df$group <- label

  return(df)
}

#Count genes in both datasets
x_counts <- count_genes(x_dir, "Ancestral X")
neo_counts <- count_genes(neo_dir, "neoX")

gene_counts <- bind_rows(x_counts, neo_counts)

print(gene_counts)

#Plot
gene_counts$group <- factor(gene_counts$group, levels = c("Ancestral X", "neoX"))
gene_counts$xpos <- ifelse(gene_counts$group == "Ancestral X", 1, 1.6)

p <- ggplot(gene_counts, aes(x = xpos, y = x_genes, group = group)) +

  geom_boxplot(width = 0.25, outlier.shape = NA) +

  geom_jitter(width = 0.05, size = 2, alpha = 0.7) +

  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.25,
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
    y = "Number of genes",
    title = "Contrast X and Neo-X genes"
  ) +

  theme_classic(base_size = 14)

ggsave("X_vs_neoX_gene_counts_boxplot.pdf", p, width = 6, height = 6)
ggsave("X_vs_neoX_gene_counts_boxplot.png", p, width = 6, height = 6, dpi = 300)
