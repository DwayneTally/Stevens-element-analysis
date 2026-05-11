#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ape)
  library(phylopairs)
  library(rstan)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3)
  stop("Usage: Rscript run_phylopairs_superfamily_updated_beta_final.R <index> <data_dir> <tree_file>")

sf_index  <- as.integer(args[1])
data_dir  <- args[2]
tree_file <- args[3]

superfamily_colors <- c(
  "Tenebrionoidea" = "#5a97cf",
  "Chrysomeloidea" = "#c1ddb5",
  "Caraboidea"     = "#dbb3c7",
  "Coccinelloidea" = "#f3cb50",
  "Elateroidea"    = "#80b569",
  "Staphylinoidea" = "#b7b7b7",
  "Scarabaeoidea"  = "#f8be81",
  "Curculionoidea" = "#eaa155",
  "Cerambycidae"   = "#c1ddb5",
  "Chrysomelidae"  = "#c1ddb5"
)

col_stevens <- "#5a97cf"

#single_model
file_sets <- list(
  list(sf = "Tenebrionoidea",
       a  = "Tenebrionoidea_pairwise_block_stats_with_tree_distance.csv",
       b  = NA,
       ref = "Tribolium_castaneum_GCF_000002335.3",
       ref_type = "within",
       single_model = TRUE),
  list(sf = "Chrysomeloidea",
       a  = "chrysomeloidea_subsampled/Chrysomeloidea_stevens_subsampled.csv",
       b  = "chrysomeloidea_subsampled/Chrysomeloidea_Gpolygoni_subsampled.csv",
       ref = "Galerucella_polygoni_GCA_905147045.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Caraboidea",
       a  = "Caraboidea_pairwis_block_stats_with_tree_distance.csv",
       b  = "Caraboidea_A_parallelepipedus_distance.csv",
       ref = "Abax_parallelepipedus_GCA_964197645.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Coccinelloidea",
       a  = "Coccinelloidea_pairwis_block_stats_with_tree_distance.csv",
       b  = "Coccinelloidea_H_sedecimguttata_distance.csv",
       ref = "Harmonia_sedecimguttata_GCA_964030865.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Elateroidea",
       a  = "Elateroidea_pairwise_block_stats_with_tree_distance.csv",
       b  = "Elat_M_seriepunctatus_distance.csv",
       ref = "Melanotus_seriepunctatus_GCA_964228145.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Staphylinoidea",
       a  = "Staphylinoidea_pairwis_block_stats_with_tree_distance.csv",
       b  = "Staphylinoidea_L_lunulatus_distance.csv",
       ref = "Lordithon_lunulatus_GCA_963942505.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Scarabaeoidea",
       a  = "Scarabaeoidea_pairwis_block_stats_with_tree_distance.csv",
       b  = "Scarab_M_prodromus_distance.csv",
       ref = "Melinopterus_prodromus_GCA_963920655.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Curculionoidea",
       a  = "Curculionoidea_pairwis_block_stats_with_tree_distance.csv",
       b  = "Curculionidea_O_rusci_distance.csv",
       ref = "Orchestes_rusci_GCA_947577165.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Cerambycidae",
       a  = "chrysomeloidea_split_subsets/Cerambycidae_stevens.csv",
       b  = "chrysomeloidea_split_subsets/Cerambycidae_superfamily.csv",
       ref = "Galerucella_polygoni_GCA_905147045.1",
       ref_type = "within",
       single_model = FALSE),
  list(sf = "Chrysomelidae",
       a  = "chrysomelidae_subsampled_seed_52/Chrysomelidae_stevens_subsampled.csv",
       b  = "chrysomelidae_subsampled_seed_52/Chrysomelidae_superfamily_subsampled.csv",
       ref = "Galerucella_polygoni_GCA_905147045.1",
       ref_type = "within",
       single_model = FALSE)
)

if (sf_index < 1 || sf_index > length(file_sets))
  stop("superfamily_index must be between 1 and ", length(file_sets))

fs <- file_sets[[sf_index]]

out_dir <- file.path("results_sf_compare", fs$sf)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
cat("Output directory:", out_dir, "\n")

data_path <- function(f) file.path(data_dir, f)

get_xmax <- function(sf_name) {
  return(450)
}

#load calibrated tree
tree <- read.tree(tree_file)
cat("Tree tips:", length(tree$tip.label),
    "| Ultrametric:", is.ultrametric(tree), "\n")

process_dataset <- function(csv_file, tree, label, out_dir) {

  cat("\n", strrep("-", 60), "\n")
  cat(" Processing:", label, "\n")
  cat(strrep("-", 60), "\n")

  dat <- read.csv(csv_file, stringsAsFactors = FALSE)
  colnames(dat)[colnames(dat) == "genome1"] <- "sp1"
  colnames(dat)[colnames(dat) == "genome2"] <- "sp2"

  cat("Rows:", nrow(dat),
      "| Species:", length(unique(c(dat$sp1, dat$sp2))), "\n")
  cat("Synteny range:", range(dat$pct_cov_pairmean), "\n")

  sp_in_data <- unique(c(dat$sp1, dat$sp2))
  missing <- setdiff(sp_in_data, tree$tip.label)
  if (length(missing) > 0)
    stop("Species missing from tree: ", paste(missing, collapse = ", "))

  tree_pruned <- keep.tip(tree, sp_in_data)
  cat("All species found in tree\n")

  coph_mat <- cophenetic(tree_pruned)
  dat$div_time <- mapply(function(s1, s2) coph_mat[s1, s2], dat$sp1, dat$sp2)
  cat("Cophenetic distance range (Ma):", range(dat$div_time), "\n")

  N <- nrow(dat)
  y_raw <- dat$pct_cov_pairmean / 100
  y <- (y_raw * (N - 1) + 0.5) / N

  cat("Calculating CP matrix...\n")
  CP <- taxapair.vcv(
    sp.pairs = dat[, c("sp1", "sp2")],
    tree     = tree_pruned,
    model    = "sq.diff",
    spnames  = FALSE
  )
  cat("CP validity:\n")
  print(covmat.check(CP))
#Normalizing covariance matrix to have a trace of 1
  CP_scaled <- CP / sum(diag(CP))
  cat("CP trace before scaling:", sum(diag(CP)), "\n")
  cat("CP trace after  scaling:", sum(diag(CP_scaled)), "\n")

  predictor <- log1p(dat$div_time)

  cat("\nBeta mixed model (beta.mm)\n")
  fit_beta_mm <- betareg.stan(
    des     = predictor,
    y       = y,
    model   = "beta.mm",
    link    = "logit",
    covmat  = CP_scaled,
    iter    = 10000,
    chains  = 4,
    cores   = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  print(fit_beta_mm[[1]])

  cat("\nStandard beta regression (beta.reg)\n")
  fit_beta_reg <- betareg.stan(
    des     = predictor,
    y       = y,
    model   = "beta.reg",
    link    = "logit",
    iter    = 10000,
    chains  = 4,
    cores   = 4,
    control = list(adapt_delta = 0.99, max_treedepth = 15)
  )
  print(fit_beta_reg[[1]])

  save_summary_csv <- function(fit, suffix) {
    tbl <- as.data.frame(fit[[1]])
    tbl$Parameter <- rownames(tbl)
    rownames(tbl) <- NULL
    tbl <- tbl[, c("Parameter", setdiff(names(tbl), "Parameter"))]
    fname <- file.path(
      out_dir,
      paste0(gsub("[^A-Za-z0-9]", "_", label), suffix, ".csv")
    )
    write.csv(tbl, file = fname, row.names = FALSE)
    cat("Summary saved to:", fname, "\n")
  }

  save_summary_csv(fit_beta_mm,  "_beta_mm_summary")
  save_summary_csv(fit_beta_reg, "_beta_reg_summary")

  rda_name <- file.path(
    out_dir,
    paste0(gsub("[^A-Za-z0-9]", "_", label), "_fit.rda")
  )
  save(fit_beta_mm, fit_beta_reg, dat, CP, CP_scaled, file = rda_name)
  cat("Models saved to:", rda_name, "\n")

  list(
    dat          = dat,
    fit_beta_mm  = fit_beta_mm,
    fit_beta_reg = fit_beta_reg,
    mm_slope  = round(fit_beta_mm[[1]]["Coef[2]", "mean"], 4),
    mm_ci_lo  = round(fit_beta_mm[[1]]["Coef[2]", "2.5%"],  3),
    mm_ci_hi  = round(fit_beta_mm[[1]]["Coef[2]", "97.5%"], 3),
    mm_sig2   = round(fit_beta_mm[[1]]["sig2_scale[1]", "mean"], 6),
    reg_slope = round(fit_beta_reg[[1]]["Coef[2]", "mean"], 4),
    reg_ci_lo = round(fit_beta_reg[[1]]["Coef[2]", "2.5%"],  3),
    reg_ci_hi = round(fit_beta_reg[[1]]["Coef[2]", "97.5%"], 3)
  )
}

beta_ribbon <- function(fit, x_seq, n_draws = 4000) {
  b0_mean <- fit[[1]]["Coef[1]", "mean"]
  b0_sd   <- fit[[1]]["Coef[1]", "sd"]
  b1_mean <- fit[[1]]["Coef[2]", "mean"]
  b1_sd   <- fit[[1]]["Coef[2]", "sd"]

  fit_line <- plogis(b0_mean + b1_mean * log1p(x_seq)) * 100

  set.seed(42)
  b0_d <- rnorm(n_draws, b0_mean, b0_sd)
  b1_d <- rnorm(n_draws, b1_mean, b1_sd)
  pred_mat <- outer(b0_d, rep(1, length(x_seq))) +
    outer(b1_d, log1p(x_seq))
  pred_pct <- plogis(pred_mat) * 100

  list(
    fit = fit_line,
    lo  = apply(pred_pct, 2, quantile, 0.025),
    hi  = apply(pred_pct, 2, quantile, 0.975)
  )
}

plot_single_model <- function(res, col_ref, sf_name, out_dir) {
  x_max <- get_xmax(sf_name)
  x_seq <- seq(0, x_max, length.out = 400)

  rib_mm  <- beta_ribbon(res$fit_beta_mm,  x_seq)
  rib_reg <- beta_ribbon(res$fit_beta_reg, x_seq)

  plot_file <- file.path(out_dir, paste0(sf_name, "_synteny_plot_single.pdf"))
  pdf(plot_file, width = 7, height = 6)

  plot(
    NA,
    xlim = c(0, x_max), ylim = c(0, 100),
    xlab = "Cophenetic distance (Mya)",
    ylab = "Syntenic coverage (%)",
    main = sf_name
  )

  polygon(
    c(x_seq, rev(x_seq)),
    c(rib_mm$lo, rev(rib_mm$hi)),
    col = adjustcolor(col_ref, alpha.f = 0.2),
    border = NA
  )

  points(
    res$dat$div_time, res$dat$pct_cov_pairmean,
    pch = 21, cex = 0.9,
    bg  = adjustcolor(col_ref, alpha.f = 0.4),
    col = "black", lwd = 0.5
  )

  lines(x_seq, rib_mm$fit,  col = col_ref, lwd = 2, lty = 1)
  lines(x_seq, rib_reg$fit, col = col_ref, lwd = 2, lty = 2)

  legend(
    "topright",
    legend = c(paste0(sf_name, " (beta.mm)"),
               paste0(sf_name, " (beta.reg)")),
    col    = c(col_ref, col_ref),
    lty    = c(1, 2),
    lwd    = 2,
    bty    = "n",
    cex    = 0.8
  )

  dev.off()
  cat("Single-model plot saved:", plot_file, "\n")
}

plot_mm_only <- function(res_st, res_ref, col_ref, sf_name, out_dir) {
  x_max <- get_xmax(sf_name)
  x_seq <- seq(0, x_max, length.out = 400)

  rib_st  <- beta_ribbon(res_st$fit_beta_mm,  x_seq)
  rib_ref <- beta_ribbon(res_ref$fit_beta_mm, x_seq)

  plot_file <- file.path(out_dir, paste0(sf_name, "_synteny_plot_mm.pdf"))
  pdf(plot_file, width = 7, height = 6)

  plot(
    NA,
    xlim = c(0, x_max), ylim = c(0, 100),
    xlab = "Cophenetic distance (Mya)",
    ylab = "Syntenic coverage (%)",
    main = sf_name
  )

  polygon(c(x_seq, rev(x_seq)),
          c(rib_st$lo, rev(rib_st$hi)),
          col = adjustcolor(col_stevens, alpha.f = 0.2), border = NA)
  polygon(c(x_seq, rev(x_seq)),
          c(rib_ref$lo, rev(rib_ref$hi)),
          col = adjustcolor(col_ref, alpha.f = 0.2), border = NA)

  points(res_st$dat$div_time, res_st$dat$pct_cov_pairmean,
         pch = 21, cex = 0.9,
         bg  = adjustcolor(col_stevens, alpha.f = 0.4),
         col = "black", lwd = 0.5)
  points(res_ref$dat$div_time, res_ref$dat$pct_cov_pairmean,
         pch = 21, cex = 0.9,
         bg  = adjustcolor(col_ref, alpha.f = 0.4),
         col = "black", lwd = 0.5)

  lines(x_seq, rib_st$fit,  col = col_stevens, lwd = 2)
  lines(x_seq, rib_ref$fit, col = col_ref, lwd = 2)

  legend("topright",
         legend = c("Stevens elements", sf_name),
         pt.bg  = c(col_stevens, col_ref),
         col    = "black",
         pch    = 21, lty = 1, lwd = 2,
         bty    = "n", cex = 0.8)

  dev.off()
  cat("Plot 1 saved:", plot_file, "\n")
}

plot_mm_vs_reg <- function(res_st, res_ref, col_ref, sf_name, out_dir) {
  x_max <- get_xmax(sf_name)
  x_seq <- seq(0, x_max, length.out = 400)

  rib_st_mm   <- beta_ribbon(res_st$fit_beta_mm,   x_seq)
  rib_ref_mm  <- beta_ribbon(res_ref$fit_beta_mm,  x_seq)
  rib_st_reg  <- beta_ribbon(res_st$fit_beta_reg,  x_seq)
  rib_ref_reg <- beta_ribbon(res_ref$fit_beta_reg, x_seq)

  plot_file <- file.path(out_dir, paste0(sf_name, "_synteny_plot_compare.pdf"))
  pdf(plot_file, width = 7, height = 6)

  plot(
    NA,
    xlim = c(0, x_max), ylim = c(0, 100),
    xlab = "Cophenetic distance (Mya)",
    ylab = "Syntenic coverage (%)",
    main = sf_name
  )

  polygon(c(x_seq, rev(x_seq)),
          c(rib_st_mm$lo, rev(rib_st_mm$hi)),
          col = adjustcolor(col_stevens, alpha.f = 0.2), border = NA)
  polygon(c(x_seq, rev(x_seq)),
          c(rib_ref_mm$lo, rev(rib_ref_mm$hi)),
          col = adjustcolor(col_ref, alpha.f = 0.2), border = NA)

  points(res_st$dat$div_time, res_st$dat$pct_cov_pairmean,
         pch = 21, cex = 0.9,
         bg  = adjustcolor(col_stevens, alpha.f = 0.4),
         col = "black", lwd = 0.5)
  points(res_ref$dat$div_time, res_ref$dat$pct_cov_pairmean,
         pch = 21, cex = 0.9,
         bg  = adjustcolor(col_ref, alpha.f = 0.4),
         col = "black", lwd = 0.5)

  lines(x_seq, rib_st_mm$fit,  col = col_stevens, lwd = 2, lty = 1)
  lines(x_seq, rib_ref_mm$fit, col = col_ref,     lwd = 2, lty = 1)

  lines(x_seq, rib_st_reg$fit,  col = col_stevens, lwd = 2, lty = 2)
  lines(x_seq, rib_ref_reg$fit, col = col_ref,     lwd = 2, lty = 2)

  legend("topright",
         legend = c("Stevens",
                    paste0(sf_name),
                    "Stevens (beta.reg)",
                    paste0(sf_name, " (beta.reg)")),
         col    = c(col_stevens, col_ref, col_stevens, col_ref),
         lty    = c(1, 1, 2, 2),
         lwd    = 2,
         bty    = "n", cex = 0.8)

  dev.off()
  cat("Plot 2 saved:", plot_file, "\n")
}

col_ref   <- superfamily_colors[fs$sf]
ref_label <- sub("_GC[AF]_.*$", "", fs$ref)

if (isTRUE(fs$single_model)) {

  label_a <- paste0(fs$sf, " — within-superfamily")

  res_single <- process_dataset(data_path(fs$a), tree, label_a, out_dir)

  plot_single_model(res_single, col_ref, fs$sf, out_dir)

  combined_rda <- file.path(out_dir, paste0(fs$sf, "_results.rda"))
  save(res_single, file = combined_rda)
  cat("Combined results saved to:", combined_rda, "\n")

  summary_file <- file.path(out_dir, paste0(fs$sf, "_summary.txt"))
  sink(summary_file)
  cat("Superfamily:", fs$sf, "\n\n")
  cat("Within-superfamily model\n")
  cat("  --- beta.mm ---\n")
  cat("  slope:      ", res_single$mm_slope, "\n")
  cat("  95% CI:     [", res_single$mm_ci_lo, ",", res_single$mm_ci_hi, "]\n")
  cat("  sig2_scale: ", res_single$mm_sig2, "\n")
  cat("  --- beta.reg ---\n")
  cat("  slope:      ", res_single$reg_slope, "\n")
  cat("  95% CI:     [", res_single$reg_ci_lo, ",", res_single$reg_ci_hi, "]\n")
  sink()
  cat("Summary written to:", summary_file, "\n")

} else {

  label_a <- paste0(fs$sf, " — T. castaneum reference")
  label_b <- paste0(fs$sf, " — ", ref_label, " reference")

  res_stevens <- process_dataset(data_path(fs$a), tree, label_a, out_dir)
  res_ref     <- process_dataset(data_path(fs$b), tree, label_b, out_dir)

  plot_mm_only(res_stevens, res_ref, col_ref, fs$sf, out_dir)
  plot_mm_vs_reg(res_stevens, res_ref, col_ref, fs$sf, out_dir)

  combined_rda <- file.path(out_dir, paste0(fs$sf, "_results.rda"))
  save(res_stevens, res_ref, file = combined_rda)
  cat("Combined results saved to:", combined_rda, "\n")

  summary_file <- file.path(out_dir, paste0(fs$sf, "_summary.txt"))
  sink(summary_file)
  cat("Superfamily:", fs$sf, "\n\n")

  cat("T. castaneum ref\n")
  cat("  --- beta.mm ---\n")
  cat("  slope:      ", res_stevens$mm_slope, "\n")
  cat("  95% CI:     [", res_stevens$mm_ci_lo, ",", res_stevens$mm_ci_hi, "]\n")
  cat("  sig2_scale: ", res_stevens$mm_sig2, "\n")
  cat("  --- beta.reg ---\n")
  cat("  slope:      ", res_stevens$reg_slope, "\n")
  cat("  95% CI:     [", res_stevens$reg_ci_lo, ",", res_stevens$reg_ci_hi, "]\n\n")

  cat(ref_label, "ref\n")
  cat("  --- beta.mm ---\n")
  cat("  slope:      ", res_ref$mm_slope, "\n")
  cat("  95% CI:     [", res_ref$mm_ci_lo, ",", res_ref$mm_ci_hi, "]\n")
  cat("  sig2_scale: ", res_ref$mm_sig2, "\n")
  cat("  --- beta.reg ---\n")
  cat("  slope:      ", res_ref$reg_slope, "\n")
  cat("  95% CI:     [", res_ref$reg_ci_lo, ",", res_ref$reg_ci_hi, "]\n")
  sink()
  cat("Summary written to:", summary_file, "\n")
}

cat("\nDone:", fs$sf, "\n")
