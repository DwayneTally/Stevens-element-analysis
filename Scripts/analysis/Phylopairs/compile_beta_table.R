#!/usr/bin/env Rscript
#Reads all _beta_mm_summary.csv and _beta_ols_summary.csv files

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1)
  stop("Usage: Rscript compile_clean_table.R <results_root> [output_csv]")

results_root <- args[1]
out_csv      <- if (length(args) >= 2) args[2] else
                file.path(results_root, "beta_table_clean.csv")

sf_root <- file.path(results_root, "results_sf_compare")

superfamilies <- c(
  "Tenebrionoidea", "Chrysomeloidea", "Caraboidea",
  "Coccinelloidea", "Elateroidea",    "Staphylinoidea",
  "Scarabaeoidea",  "Curculionoidea", "Cerambycidae",
  "Chrysomelidae"
)

extract_params <- function(csv_path) {
  if (!file.exists(csv_path)) return(NULL)

  tbl <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = TRUE)
  rownames(tbl) <- tbl$Parameter

  get_val <- function(param, col_variants) {
    if (!param %in% rownames(tbl)) return(NA)
    for (col in col_variants)
      if (col %in% colnames(tbl))
        return(tbl[param, col])
    return(NA)
  }

  mean_col <- c("mean")
  lo_col   <- c("X2.5.", "2.5.", "X2.5%")
  hi_col   <- c("X97.5.", "97.5.", "X97.5%")
  rhat_col <- c("Rhat", "rhat")
  neff_col <- c("n_eff", "neff", "Bulk_ESS")

  phi_param <- rownames(tbl)[grepl("^phi", rownames(tbl))][1]
  if (is.na(phi_param)) phi_param <- "phi"

  # Back-transform intercept for readability: plogis(logit) * 100
  intercept_logit <- get_val("Coef[1]", mean_col)
  intercept_pct   <- if (!is.na(intercept_logit))
                       round(plogis(intercept_logit) * 100, 2)
                     else NA

  data.frame(
    intercept_mean  = intercept_pct,
    intercept_Rhat  = get_val("Coef[1]", rhat_col),
    slope_mean      = get_val("Coef[2]", mean_col),
    slope_lo95      = get_val("Coef[2]", lo_col),
    slope_hi95      = get_val("Coef[2]", hi_col),
    slope_Rhat      = get_val("Coef[2]", rhat_col),
    slope_neff      = get_val("Coef[2]", neff_col),
    sig2_scale_mean = get_val("sig2_scale[1]", mean_col),
    sig2_scale_Rhat = get_val("sig2_scale[1]", rhat_col),
    phi_mean        = get_val(phi_param, mean_col),
    stringsAsFactors = FALSE
  )
}

parse_label <- function(fname) {
  base    <- tools::file_path_sans_ext(basename(fname))
  model   <- if (grepl("beta_mm",  base)) "beta.mm"  else "beta.reg"
  dataset <- if (grepl("castaneum|T__castaneum", base, ignore.case = TRUE))
               "Stevens (T. castaneum ref)"
             else
               "Within-SF ref"
  list(model = model, dataset = dataset)
}

all_rows <- list()

for (sf in superfamilies) {
  sf_dir <- file.path(sf_root, sf)
  if (!dir.exists(sf_dir)) {
    cat("Skipping (not found):", sf, "\n"); next
  }

  csvs <- list.files(sf_dir,
                     pattern = "_(beta_mm|beta_ols)_summary\\.csv$",
                     full.names = TRUE)
  if (length(csvs) == 0) {
    cat("No CSVs found for:", sf, "\n"); next
  }

  cat(sf, "—", length(csvs), "CSVs\n")

  parsed_list <- lapply(csvs, parse_label)
  keys        <- sapply(parsed_list,
                        function(p) paste(p$dataset, p$model, sep = "|"))
  mtimes      <- file.info(csvs)$mtime

  deduped <- c()
  for (key in unique(keys)) {
    idx  <- which(keys == key)
    pick <- idx[which.max(mtimes[idx])]
    if (length(idx) > 1)
      cat("  Deduped:", key, "— kept", basename(csvs[pick]), "\n")
    deduped <- c(deduped, csvs[pick])
  }

  for (csv_path in deduped) {
    parsed <- parse_label(csv_path)
    params <- extract_params(csv_path)
    if (is.null(params)) next

    # Get n_pairs and n_species from rda if available
    n_pairs <- NA; n_species <- NA
    rda <- file.path(sf_dir, paste0(sf, "_results.rda"))
    if (file.exists(rda)) {
      env <- new.env()
      tryCatch({
        load(rda, envir = env)
        res_obj <- if (grepl("castaneum", parsed$dataset, ignore.case = TRUE)) {
          if (exists("res_stevens", envir = env)) env$res_stevens else NULL
        } else {
          if (exists("res_ref", envir = env)) env$res_ref else NULL
        }
        if (!is.null(res_obj) && !is.null(res_obj$dat)) {
          n_pairs   <- nrow(res_obj$dat)
          n_species <- length(unique(c(res_obj$dat$sp1, res_obj$dat$sp2)))
        }
      }, error = function(e) NULL)
    }

    row <- data.frame(
      superfamily = sf,
      dataset     = parsed$dataset,
      model       = parsed$model,
      n_pairs     = n_pairs,
      n_species   = n_species,
      stringsAsFactors = FALSE
    )
    row <- cbind(row, params)
    all_rows[[length(all_rows) + 1]] <- row
    cat("  ", parsed$dataset, "|", parsed$model, "\n")
  }
}

if (length(all_rows) == 0) {
  cat("\nNo data found.\n"); quit(status = 1)
}

master <- do.call(rbind, all_rows)
rownames(master) <- NULL

num_cols <- sapply(master, is.numeric)
master[num_cols] <- lapply(master[num_cols], round, 4)

#Order by superfamily, then dataset, then model
master <- master[order(master$superfamily,
                       master$dataset,
                       master$model), ]

write.csv(master, file = out_csv, row.names = FALSE)

#Preview
print(master[, c("superfamily", "dataset", "model",
                  "intercept_mean", "intercept_Rhat",
                  "slope_mean", "slope_lo95", "slope_hi95",
                  "sig2_scale_mean")],
      row.names = FALSE)
