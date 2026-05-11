library(tidyverse)

bed_dir <- "Genes_on_X_noNeoX/bed"

bed_files <- list.files(bed_dir, pattern = "\\.bed$", full.names = TRUE)

if (length(bed_files) == 0) {
  stop("No .bed files found in: ", bed_dir)
}

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
)

gene_counts$genome <- sub("\\.bed$", "", gene_counts$genome)
gene_counts$group <- "X"

print(gene_counts)

p <- ggplot(gene_counts, aes(x = group, y = x_genes)) +
  geom_boxplot(width = 0.25, outlier.shape = NA) +
  geom_jitter(width = 0.08, size = 2, alpha = 0.7) +
  labs(
    x = "",
    y = "Genes on X chromosome",
    title = "X-linked genes"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

#print(p)

ggsave("X_gene_counts_boxplot.pdf", p, width = 5, height = 6)
ggsave("X_gene_counts_boxplot.png", p, width = 5, height = 6, dpi = 300)
