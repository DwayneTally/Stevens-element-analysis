library(data.table)

stevens_file <- "Chrysomeloidea_pairwis_block_stats_with_tree_distance.csv"
superfamily_file <- "Chrysomeliodea_G_Polygoni_distance.csv"

out_dir <- "chrysomeloidea_split_subsets"
dir.create(out_dir, showWarnings = FALSE)


cerambycidae_raw <- c(
  "Pogonocherus_hispidulus_Cerambycidae_GCA_963924545.1",
  "Monochamus_alternatus_Cerambycidae_GCA_037114965.1",
  "Monochamus_saltuarius_Cerambycidae_GCA_025584915.1",
  "Arhopalus_rusticus_Cerambycidae_GCA_964213965.1",
  "Tetropium_fuscum_Cerambycidae_GCA_964058775.1",
  "Rhagium_mordax_Cerambycidae_GCA_963680705.1",
  "Stenurella_melanura_Cerambycidae_GCA_963583905.1",
  "Rutpela_maculata_Cerambycidae_GCA_936432065.2",
  "Leptura_quadrifasciata_Cerambycidae_GCA_963675555.1",
  "Stictoleptura_scutellata_Cerambycidae_GCA_964212005.1"
)

cerambycidae <- sub("_Cerambycidae_", "_", cerambycidae_raw)

cat("Cerambycidae taxa:\n")
print(cerambycidae)

stev <- fread(stevens_file)
sup  <- fread(superfamily_file)

required_cols <- c("genome1", "genome2")
if (!all(required_cols %in% names(stev))) {
  stop("Stevens file missing genome1/genome2")
}
if (!all(required_cols %in% names(sup))) {
  stop("Superfamily file missing genome1/genome2")
}

all_stev_taxa <- sort(unique(c(stev$genome1, stev$genome2)))
all_sup_taxa  <- sort(unique(c(sup$genome1, sup$genome2)))
all_taxa <- sort(unique(c(all_stev_taxa, all_sup_taxa)))

cat("\nNumber of taxa across both files:", length(all_taxa), "\n")

chrysomelidae <- setdiff(all_taxa, c(cerambycidae, tribolium))

cat("\nInferred Chrysomelidae taxa:\n")
print(chrysomelidae)

missing_cerambycidae <- setdiff(cerambycidae, all_taxa)
if (length(missing_cerambycidae) > 0) {
  cat("\nWARNING: These Cerambycidae taxa were not found in the CSVs:\n")
  print(missing_cerambycidae)
}

subset_stevens <- function(dt, focal_taxa, tribolium) {
  allowed <- unique(c(focal_taxa, tribolium))
  out <- dt[genome1 %in% allowed & genome2 %in% allowed]

  # remove self comparisons
  out <- out[genome1 != genome2]

  # deduplicate A-B vs B-A
  out[, g1 := pmin(genome1, genome2)]
  out[, g2 := pmax(genome1, genome2)]
  out <- unique(out, by = c("g1", "g2"))
  out[, c("g1", "g2") := NULL]

  out
}

subset_superfamily <- function(dt, focal_taxa) {
  out <- dt[genome1 %in% focal_taxa & genome2 %in% focal_taxa]

  #remove self comparisons
  out <- out[genome1 != genome2]

  out[, g1 := pmin(genome1, genome2)]
  out[, g2 := pmax(genome1, genome2)]
  out <- unique(out, by = c("g1", "g2"))
  out[, c("g1", "g2") := NULL]

  out
}
#make subsets
cerambycidae_stevens <- subset_stevens(stev, cerambycidae, tribolium)
cerambycidae_superfamily <- subset_superfamily(sup, cerambycidae)

chrysomelidae_stevens <- subset_stevens(stev, chrysomelidae, tribolium)
chrysomelidae_superfamily <- subset_superfamily(sup, chrysomelidae)

cat("\nCerambycidae stevens rows:", nrow(cerambycidae_stevens), "\n")
cat("Cerambycidae superfamily rows:", nrow(cerambycidae_superfamily), "\n")

cat("\nChrysomelidae stevens rows:", nrow(chrysomelidae_stevens), "\n")
cat("Chrysomelidae superfamily rows:", nrow(chrysomelidae_superfamily), "\n")

fwrite(cerambycidae_stevens,
       file.path(out_dir, "Cerambycidae_stevens.csv"))

fwrite(cerambycidae_superfamily,
       file.path(out_dir, "Cerambycidae_superfamily.csv"))

fwrite(chrysomelidae_stevens,
       file.path(out_dir, "Chrysomelidae_stevens.csv"))

fwrite(chrysomelidae_superfamily,
       file.path(out_dir, "Chrysomelidae_superfamily.csv"))

fwrite(data.table(genome = cerambycidae),
       file.path(out_dir, "Cerambycidae_taxa.csv"))

fwrite(data.table(genome = chrysomelidae),
       file.path(out_dir, "Chrysomelidae_taxa.csv"))

cat("\nDone.\n")
