library(data.table)
library(ape)
library(phytools)
library(phangorn)

#Cai MCMCtree output files
cai_IR1_file  <- "IR/chain1/FigTree.tre"
cai_IR2_file  <- "IR/chain2/FigTree.tre"
cai_AC1_file  <- "AC/chain1/FigTree.tre"
cai_AC2_file  <- "AC/chain2/FigTree.tre"

tree_file     <- "Tree_branch_lengths.txt"
meta_file     <- "genome_superfamily.tsv"
families_file <- "Beetle_genome_families.txt"
out_prefix    <- "beetle_timetree_Cai2022"

#MCMCtree stores branch lengths in units of 100 Mya
#Branch lengths multiplied by 100 to convert to Mya

parse_mcmctree <- function(filepath) {
  cat("Reading:", filepath, "\n")
  raw <- readLines(filepath)
  
  tree_line <- raw[grep("UTREE", raw)]
  if (length(tree_line) == 0) stop(paste("No UTREE found in:", filepath))
  
  tree_line <- gsub(".*UTREE \\d+ = ", "", tree_line)
  tree_line <- gsub("\\[&95%=\\{[^}]*\\}\\]", "", tree_line)
  
  writeLines(tree_line, "temp_cai_tree.nwk")
  tree <- read.tree("temp_cai_tree.nwk")
  
  tree$edge.length <- tree$edge.length * 100
  
  cat("  Tips:", length(tree$tip.label), "\n")
  cat("  Root age (Ma):", round(max(node.depth.edgelength(tree)), 1), "\n\n")
  return(tree)
}

tree_IR1 <- parse_mcmctree(cai_IR1_file)
tree_IR2 <- parse_mcmctree(cai_IR2_file)
tree_AC1 <- parse_mcmctree(cai_AC1_file)
tree_AC2 <- parse_mcmctree(cai_AC2_file)

print(tree_IR1$tip.label)

get_node_ages <- function(tree) {
  max(node.depth.edgelength(tree)) - node.depth.edgelength(tree)
}

ages_IR1 <- get_node_ages(tree_IR1)
ages_IR2 <- get_node_ages(tree_IR2)
ages_AC1 <- get_node_ages(tree_AC1)
ages_AC2 <- get_node_ages(tree_AC2)

#pattern: grep string matching relevant tip labels in Cai tree
#Returns a one-row data.frame with ages from all four chains,
#plus mean, min, max and convergence assessment

get_clade_ages <- function(clade_name, pattern) {
  
  tips_IR1 <- grep(pattern, tree_IR1$tip.label, value = TRUE)
  tips_IR2 <- grep(pattern, tree_IR2$tip.label, value = TRUE)
  tips_AC1 <- grep(pattern, tree_AC1$tip.label, value = TRUE)
  tips_AC2 <- grep(pattern, tree_AC2$tip.label, value = TRUE)
  
  if (length(tips_IR1) < 2) {
    cat(sprintf("%-22s — fewer than 2 tips matched (pattern: '%s')\n",
                clade_name, pattern))
    return(NULL)
  }
  
  node_IR1 <- getMRCA(tree_IR1, tips_IR1)
  node_IR2 <- getMRCA(tree_IR2, tips_IR2)
  node_AC1 <- getMRCA(tree_AC1, tips_AC1)
  node_AC2 <- getMRCA(tree_AC2, tips_AC2)
  
  age_IR1 <- round(ages_IR1[node_IR1], 1)
  age_IR2 <- round(ages_IR2[node_IR2], 1)
  age_AC1 <- round(ages_AC1[node_AC1], 1)
  age_AC2 <- round(ages_AC2[node_AC2], 1)
  
  all_ages   <- c(age_IR1, age_IR2, age_AC1, age_AC2)
  mean_age   <- round(mean(all_ages), 1)
  min_age    <- min(all_ages)
  max_age    <- max(all_ages)
  ir_ac_diff <- abs(mean(c(age_IR1, age_IR2)) - mean(c(age_AC1, age_AC2)))
  confidence <- ifelse(ir_ac_diff <= 5,  "HIGH",
                       ifelse(ir_ac_diff <= 15, "MEDIUM", "LOW"))
  
  cat(sprintf("\n%s\n", clade_name))
  cat(sprintf("  IR chain1:   %.1f Ma\n", age_IR1))
  cat(sprintf("  IR chain2:   %.1f Ma\n", age_IR2))
  cat(sprintf("  AC chain1:   %.1f Ma\n", age_AC1))
  cat(sprintf("  AC chain2:   %.1f Ma\n", age_AC2))
  cat(sprintf("  Mean:        %.1f Ma\n", mean_age))
  cat(sprintf("  Range:       %.1f - %.1f Ma\n", min_age, max_age))
  cat(sprintf("  IR vs AC:    %.1f Ma\n", ir_ac_diff))
  cat(sprintf("  Confidence:  %s\n", confidence))
  
  data.frame(
    clade      = clade_name,
    IR_chain1  = age_IR1,
    IR_chain2  = age_IR2,
    AC_chain1  = age_AC1,
    AC_chain2  = age_AC2,
    mean_Ma    = mean_age,
    age_min    = min_age,
    age_max    = max_age,
    IR_vs_AC   = round(ir_ac_diff, 1),
    confidence = confidence,
    stringsAsFactors = FALSE
  )
}

chain_results <- list(
  get_clade_ages("Carabidae",      "Carabidae"),
  get_clade_ages("Staphylinidae",  "Staphylinidae"),
  get_clade_ages("Scarabaeidae",   "Scarabaeidae"),
  get_clade_ages("Cerambycidae",   "Cerambycidae"),
  get_clade_ages("Chrysomelidae",  "Chrysomelidae"),
  get_clade_ages("Curculionidae",  "Curculionidae"),
  get_clade_ages("Coccinellidae",  "Coccinellidae"),
  get_clade_ages("Tenebrionidae",  "Tenebrionidae"),
  get_clade_ages("Elateridae",     "Elateridae"),
  get_clade_ages("Cantharidae",    "Cantharidae"),
  get_clade_ages("Scarabaeoidea",  "Scarabaeo|Geotrup|Lucan"),
  get_clade_ages("Staphylinoidea", "Staphylin|Silph|Histeri|Leiod"),
  get_clade_ages("Elateroidea",    "Elaterid|Cantharid|Lampyrid"),
  get_clade_ages("Curculionoidea", "Curculion|Anthribid|Attelab|Apionid"),
  get_clade_ages("Cleroidea",      "Clerid|Melyrid"),
  get_clade_ages("Tenebrionoidea", "Tenebrion|Meloid|Pyrochro|Salpingid|Scraptii|Zopheri")
)

chain_df <- do.call(rbind, chain_results[!sapply(chain_results, is.null)])

print(chain_df)

write.csv(chain_df, paste0(out_prefix, "_cai_chain_ages.csv"),
          row.names = FALSE)
cat("\nChain ages saved to:",
    paste0(out_prefix, "_cai_chain_ages.csv"), "\n")

#Clean up temp
if (file.exists("temp_cai_tree.nwk")) file.remove("temp_cai_tree.nwk")

chain_age_min <- function(clade) {
  row <- chain_df[chain_df$clade == clade, ]
  if (nrow(row) == 0) stop(paste("Clade not found in chain results:", clade))
  row$age_min
}

chain_age_max <- function(clade) {
  row <- chain_df[chain_df$clade == clade, ]
  if (nrow(row) == 0) stop(paste("Clade not found in chain results:", clade))
  row$age_max
}

#Deep nodes from Cai Table 1
#Family/superfamily nodes use ages extracted from chains above

calibration_info <- data.table(
  superfamily = c(
    "Root_Strepsiptera_Coleoptera",
    "Crown_Coleoptera",
    "Adephaga",
    "Polyphaga",
    "Elateriformia",
    "Staphyliniformia",
    "Cucujiformia",
    "Elateroidea",
    "Cleroidea",
    "Curculionoidea",
    "Carabidae",
    "Scarabaeidae",
    "Coccinellidae",
    "Chrysomelidae",
    "Cerambycidae",
    "Cantharidae",
    "Curculionidae",
    "Tenebrionidae"
  ),
  age_min = c(
    #Root
    350,
    #Table 1
    306, 259, 286, 246, 238, 220,
    # From chain extraction
    chain_age_min("Elateroidea"),
    chain_age_min("Cleroidea"),
    chain_age_min("Curculionoidea"),
    chain_age_min("Carabidae"),
    chain_age_min("Scarabaeidae"),
    chain_age_min("Coccinellidae"),
    chain_age_min("Chrysomelidae"),
    chain_age_min("Cerambycidae"),
    chain_age_min("Cantharidae"),
    chain_age_min("Curculionidae"),
    chain_age_min("Tenebrionidae")
  ),
  age_max = c(
    #Root
    400,
    #Table 1
    322, 288, 307, 269, 258, 249,
    # From chain extraction
    chain_age_max("Elateroidea"),
    chain_age_max("Cleroidea"),
    chain_age_max("Curculionoidea"),
    chain_age_max("Carabidae"),
    chain_age_max("Scarabaeidae"),
    chain_age_max("Coccinellidae"),
    chain_age_max("Chrysomelidae"),
    chain_age_max("Cerambycidae"),
    chain_age_max("Cantharidae"),
    chain_age_max("Curculionidae"),
    chain_age_max("Tenebrionidae")
  )
)

cat("Calibration table built:\n")
print(calibration_info)

tr <- read.tree(tree_file)
if (is.null(tr)) stop("Could not read tree: ", tree_file)

tr$tip.label <- gsub("^busco_", "", tr$tip.label)
tr$tip.label <- gsub("_genomic$", "", tr$tip.label)

cat("Loaded tree with", length(tr$tip.label), "tips\n")
cat("Is binary:   ", is.binary(tr), "\n")
cat("Ultrametric: ", is.ultrametric(tr), "\n")
cat("Outgroup:    ", tr$tip.label[1], "\n\n")

if (any(tr$edge.length <= 0)) {
  n_fixed <- sum(tr$edge.length <= 0)
  cat("Fixing", n_fixed, "zero/negative branch lengths -> 1e-6\n")
  tr$edge.length[tr$edge.length <= 0] <- 1e-6
}

#load TSV file
meta <- fread(meta_file)
setnames(meta, trimws(names(meta)))

if (!all(c("genome", "superfamily") %in% names(meta))) {
  stop("Metadata file must have columns: genome and superfamily")
}

meta[, genome      := trimws(sub("^\ufeff", "", genome))]
meta[, superfamily := trimws(sub("^\ufeff", "", superfamily))]
meta <- unique(meta[, .(genome, superfamily)])

cat("Loaded metadata for", nrow(meta), "genomes\n")
cat("Groups in TSV:\n")
print(sort(unique(meta$superfamily)))

#Parse family-level metadata from Beetle_genome_families.txt
raw_lines <- readLines(families_file)
raw_lines <- trimws(raw_lines[nchar(trimws(raw_lines)) > 0])

parse_family_entry <- function(x) {
  parts   <- strsplit(x, "_")[[1]]
  gca_pos <- which(parts %in% c("GCA", "GCF"))
  if (length(gca_pos) == 0) return(NULL)
  gca_pos   <- gca_pos[1]
  family    <- parts[gca_pos - 1]
  genus_sp  <- paste(parts[1:(gca_pos - 2)], collapse = "_")
  accession <- paste(parts[gca_pos:length(parts)], collapse = "_")
  genome_id <- paste0(genus_sp, "_", accession)
  data.table(genome = genome_id, superfamily = family)
}

family_meta <- rbindlist(lapply(raw_lines, parse_family_entry), fill = TRUE)
family_meta <- family_meta[!is.na(genome)]

cat("\nParsed", nrow(family_meta), "entries from", families_file, "\n")
cat("Families found:\n")
print(sort(unique(family_meta$superfamily)))

#Root:Xenos + Liopterus MRCA = root node
root_group <- data.table(
  genome = c(
    grep("Xenos",     tr$tip.label, value = TRUE),
    grep("Liopterus", tr$tip.label, value = TRUE)[1]
  ),
  superfamily = "Root_Strepsiptera_Coleoptera"
)

#Crown Coleoptera: Liopterus + Tenebrio
coleoptera_group <- data.table(
  genome = c(
    grep("Liopterus",      tr$tip.label, value = TRUE)[1],
    grep("Tenebrio_molitor", tr$tip.label, value = TRUE)[1]
  ),
  superfamily = "Crown_Coleoptera"
)

#Adephaga: Carabidae + Dytiscidae
adephaga_group <- rbind(
  meta[superfamily == "Caraboidea",
       .(genome, superfamily)][, superfamily := "Adephaga"],
  family_meta[superfamily == "Dytiscidae",
              .(genome, superfamily)][, superfamily := "Adephaga"]
)

#Polyphaga: all non-Adephaga, no Strepsiptera
polyphaga_group <- meta[
  !superfamily %in% c("Caraboidea"),
  .(genome, superfamily)
][, superfamily := "Polyphaga"]
polyphaga_group <- polyphaga_group[
  !grepl("Xenos|Liopterus", genome, ignore.case = TRUE)]

#Elateriformia: Elateroidea + Buprestidae + Dascillidae + Elmidae
elateriformia_group <- rbind(
  meta[superfamily == "Elateroidea",
       .(genome, superfamily)][, superfamily := "Elateriformia"],
  family_meta[superfamily %in% c("Buprestidae", "Dascillidae", "Elmidae"),
              .(genome, superfamily)][, superfamily := "Elateriformia"]
)

#Staphyliniformia: Staphylinoidea + Scarabaeoidea
staphyliniformia_group <- meta[
  superfamily %in% c("Staphylinoidea", "Scarabaeoidea"),
  .(genome, superfamily)
][, superfamily := "Staphyliniformia"]

#Cucujiformia: Tenebrionoidea + Cleroidea + Coccinelloidea + chrysomeloidea + Curculionoidea
cucujiformia_group <- meta[
  superfamily %in% c("Tenebrionidae", "Cleroidea",
                     "Coccinelloidea", "chrysomeloidea",
                     "Curculionoidea"),
  .(genome, superfamily)
][, superfamily := "Cucujiformia"]

tenebrionidae_group <- family_meta[
  superfamily == "Tenebrionidae",
  .(genome, superfamily)
]

carabidae_group <- family_meta[
  superfamily == "Carabidae",
  .(genome, superfamily)
]

#Other family-level groups from families file
family_groups <- rbind(
  family_meta[superfamily %in% c(
    "Coccinellidae", "Chrysomelidae", "Cerambycidae",
    "Cantharidae",  "Scarabaeidae",  "Curculionidae"
  )],
  tenebrionidae_group,
  carabidae_group,
  root_group,
  coleoptera_group,
  adephaga_group,
  polyphaga_group,
  elateriformia_group,
  staphyliniformia_group,
  cucujiformia_group
)

#Combine TSV groups + all extra groups
meta_combined <- rbind(meta, family_groups)
meta_combined <- meta_combined[genome %in% tr$tip.label]
print(meta_combined[, .N, by = superfamily][order(superfamily)])

calibration_rows <- list()

for (i in seq_len(nrow(calibration_info))) {
  
  sf      <- calibration_info$superfamily[i]
  age_min <- calibration_info$age_min[i]
  age_max <- calibration_info$age_max[i]
  
  sf_tips_all  <- meta_combined[superfamily == sf, genome]
  sf_tips      <- unique(intersect(sf_tips_all, tr$tip.label))
  missing_tips <- setdiff(sf_tips_all, tr$tip.label)
  
  cat("Group:", sf, "\n")
  cat("Metadata tips:", length(sf_tips_all), "\n")
  cat("Tips found in tree:", length(sf_tips), "\n")
  
  if (length(missing_tips) > 0) {
    cat("Missing from tree:", length(missing_tips), "\n")
    print(missing_tips)
  }
  
  if (length(sf_tips) < 2) {
    cat("Skipping", sf, "- too few tips found in tree\n")
    next
  }
  
  node <- findMRCA(tr, sf_tips)
  
  if (is.null(node) || is.na(node)) {
    cat("Skipping", sf, "- could not find MRCA\n")
    next
  }
  
  clade_tips <- extract.clade(tr, node)$tip.label
  extra_tips <- setdiff(clade_tips, sf_tips)
  
  cat("MRCA node:", node, "\n")
  cat("Tips in MRCA clade:", length(clade_tips), "\n")
  cat("Extra tips inside clade:", length(extra_tips), "\n")
  
  if (length(extra_tips) > 0) {
    cat("WARNING:", sf, "may not be monophyletic\n")
    print(head(extra_tips, 10))
  }
  
  calibration_rows[[length(calibration_rows) + 1]] <- data.table(
    superfamily           = sf,
    node                  = node,
    age_min               = age_min,
    age_max               = age_max,
    source                = "Cai et al. 2022",
    n_tips_metadata       = length(sf_tips_all),
    n_tips_used           = length(sf_tips),
    n_tips_in_clade       = length(clade_tips),
    n_extra_tips_in_clade = length(extra_tips)
  )
}

calibration_nodes <- rbindlist(calibration_rows, fill = TRUE)

if (nrow(calibration_nodes) == 0) stop("No usable calibration nodes found.")

#Remove duplicate nodes, keep first (deeper) calibration
if (any(duplicated(calibration_nodes$node))) {
  dups <- calibration_nodes[
    duplicated(node) | duplicated(node, fromLast = TRUE)]
  print(dups[, .(superfamily, node, age_min, age_max)])
  calibration_nodes <- calibration_nodes[!duplicated(node)]
}

cat("Final calibration nodes:\n")
print(calibration_nodes[, .(superfamily, node, age_min, age_max, source)])

fwrite(calibration_nodes, paste0(out_prefix, "_calibration_nodes.csv"))
cat("✓ Calibration nodes saved to:",
    paste0(out_prefix, "_calibration_nodes.csv"), "\n")

calibration <- makeChronosCalib(
  tr,
  node        = calibration_nodes$node,
  age.min     = calibration_nodes$age_min,
  age.max     = calibration_nodes$age_max,
  soft.bounds = TRUE
)

cat("Chronos calibration object:\n")
print(calibration)

cat("\nRunning chronos on full tree...\n")

tr_time <- chronos(
  tr,
  calibration = calibration,
  quiet       = FALSE
)

cat("\nUltrametric after chronos:", is.ultrametric(tr_time), "\n")

write.tree(tr_time, file = paste0(out_prefix, ".newick"))
cat("✓ Timetree saved to:", paste0(out_prefix, ".newick"), "\n")

#Recompute node ages on ultrametric tree
depths   <- node.depth.edgelength(tr_time)
root_age <- max(depths)
cat("Root age (Ma):", round(root_age, 1), "\n")

ultra_rows <- list()
for (i in seq_len(nrow(calibration_info))) {
  sf      <- calibration_info$superfamily[i]
  sf_tips <- unique(intersect(
    meta_combined[superfamily == sf, genome],
    tr_time$tip.label
  ))
  if (length(sf_tips) < 2) next
  node <- findMRCA(tr_time, sf_tips)
  if (is.null(node) || is.na(node)) next
  ultra_rows[[length(ultra_rows) + 1]] <- data.table(
    superfamily = sf,
    node        = node,
    age_min     = calibration_info$age_min[i],
    age_max     = calibration_info$age_max[i],
    node_age_ma = round(root_age - depths[node], 2),
    n_tips_used = length(sf_tips)
  )
}

ultra_nodes <- rbindlist(ultra_rows, fill = TRUE)

cat("\nEstimated ages of calibrated nodes\n")
cat(sprintf("%-35s  %8s  %14s  %s\n",
            "Group", "Age (Ma)", "Target range", "Status"))
cat(strrep("-", 75), "\n")
for (i in seq_len(nrow(ultra_nodes))) {
  in_range <- ultra_nodes$node_age_ma[i] >= ultra_nodes$age_min[i] &
    ultra_nodes$node_age_ma[i] <= ultra_nodes$age_max[i]
  flag <- ifelse(in_range, "OK", "⚠ outside range")
  cat(sprintf("%-35s  %8.1f  %5.1f-%-7.1f  %s\n",
              ultra_nodes$superfamily[i],
              ultra_nodes$node_age_ma[i],
              ultra_nodes$age_min[i],
              ultra_nodes$age_max[i],
              flag))
}

fwrite(ultra_nodes, paste0(out_prefix, "_calibration_node_ages.csv"))
cat("\nNode ages saved to:",
    paste0(out_prefix, "_calibration_node_ages.csv"), "\n")

tr_time_plot <- ladderize(tr_time)
node_cols    <- seq_len(nrow(ultra_nodes)) + 1

#Annotated plot with calibration nodes marked
pdf(paste0(out_prefix, "_ultrametric_tree_with_nodes.pdf"),
    width = 12, height = 16)
plot(tr_time_plot, cex = 0.5, no.margin = TRUE,
     main = "Time-Calibrated Tree")
axisPhylo()
add.scale.bar()
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
points(x = obj$xx[ultra_nodes$node], y = obj$yy[ultra_nodes$node],
       pch = 21, bg = node_cols, cex = 2)
text(x      = obj$xx[ultra_nodes$node],
     y      = obj$yy[ultra_nodes$node],
     labels = paste0(ultra_nodes$superfamily,
                     " (target=", ultra_nodes$age_min,
                     "-", ultra_nodes$age_max,
                     " Ma; est=", round(ultra_nodes$node_age_ma, 1), " Ma)"),
     pos = 4, cex = 0.55, offset = 0.4)
legend("topleft",
       legend = paste0(ultra_nodes$superfamily,
                       " (", ultra_nodes$age_min,
                       "-", ultra_nodes$age_max, " Ma)"),
       pch = 21, pt.bg = node_cols, pt.cex = 2, bty = "n", cex = 0.6)
dev.off()
cat("Annotated plot saved to:",
    paste0(out_prefix, "_ultrametric_tree.pdf"), "\n")


cat("Output files:\n")
cat(" ", paste0(out_prefix, "_cai_chain_ages.csv"),
    "  — chain extraction summary\n")
cat(" ", paste0(out_prefix, "_calibration_nodes.csv"),
    "  — MRCA nodes used\n")
cat(" ", paste0(out_prefix, "_calibration_node_ages.csv"),
    "  — final node age estimates\n")
cat(" ", paste0(out_prefix, ".newick"),
    "  — ultrametric timetree\n")
cat(" ", paste0(out_prefix, "_ultrametric_tree_with_nodes.pdf"),
    "  — annotated plot\n")
