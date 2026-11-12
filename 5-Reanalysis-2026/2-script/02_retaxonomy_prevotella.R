################################################################################
##  02_retaxonomy_prevotella.R   (TEMPLATE - drafted Oct 2026, not yet run)
##  Re-assign the 2022 ASVs with current references so that the Prevotella split
##  (Hitch et al. 2022: Segatella / Hoylesella / Leyella / Palleniella / ...) is used
################################################################################
##  Why
##   - SILVA 138 (2021) lumps Segatella copri (ex P. copri), Leyella stercorea
##     (ex P. stercorea) and oral Prevotella into one "Prevotella" genus
##   - the 2022 Prevotella enterotype therefore pools genera with different biology
##   - the provisional crosswalk (1-input/prevotella_reassignment.csv, reassign.py)
##     aligned ASVs to 56 type strains only; this script does it properly with DADA2
##  Strategy
##   - assign every ASV against (A) SILVA 138.2 and (B) GTDB r220 (dada2-formatted)
##   - build a crosswalk 2022 genus -> new genus, weighted by reads
##   - rebuild genus tables (new / old) and re-derive enterotypes (DMM)
##  Requirements
##   - BiocManager::install(c("dada2", "DECIPHER", "DirichletMultinomial"))
##   - reference FASTA files (Zenodo) in 1-Input/db/ :
##       silva_nr99_v138.2_toGenus_trainset.fa.gz, silva_v138.2_assignSpecies.fa.gz
##       GTDB_bac120_arc53_ssu_r220_genus.fa.gz, GTDB_bac120_arc53_ssu_r220_species.fa.gz
##   - run from the repository root; ~1 h with 8 threads
################################################################################

rm(list = ls())
suppressPackageStartupMessages({ library(dada2); library(Biostrings); library(tidyverse) })
if (!dir.exists("4-Script")) stop("run from the repository root")
out_dir <- "5-Reanalysis-2026/4-output"; dir.create(out_dir, showWarnings = FALSE)

db_silva    <- "1-Input/db/silva_nr99_v138.2_toGenus_trainset.fa.gz"
db_silva_sp <- "1-Input/db/silva_v138.2_assignSpecies.fa.gz"
db_gtdb     <- "1-Input/db/GTDB_bac120_arc53_ssu_r220_genus.fa.gz"
db_gtdb_sp  <- "1-Input/db/GTDB_bac120_arc53_ssu_r220_species.fa.gz"
stopifnot(file.exists(db_silva, db_silva_sp, db_gtdb, db_gtdb_sp))

## 0. ASV sequences and counts (2022 inputs) -----------------------------------------
seqs   <- readDNAStringSet("1-Input/asv_seq.fasta")
seqtab <- read.csv("1-Input/CRC_prognosis_ASV_table.csv", row.names = 1, check.names = FALSE) %>% as.matrix() %>% t()
##  rows = samples (339), columns = ASV ids; keep ASVs with a sequence
seqtab <- seqtab[, intersect(colnames(seqtab), names(seqs))]
seq_chr <- as.character(seqs[colnames(seqtab)])
tax_old <- read.csv("1-Input/CRC_prognosis_Taxonomy.csv") %>%
  mutate(Genus = str_extract(Taxon, "g__[^;]+") %>% str_remove("g__") %>% str_trim())

## 1. Re-assignment with both databases ------------------------------------------------
set.seed(100)
tax_silva <- assignTaxonomy(seq_chr, db_silva, multithread = TRUE, minBoot = 80)
tax_silva <- addSpecies(tax_silva, db_silva_sp, allowMultiple = TRUE)
tax_gtdb  <- assignTaxonomy(seq_chr, db_gtdb,  multithread = TRUE, minBoot = 80)
tax_gtdb  <- addSpecies(tax_gtdb,  db_gtdb_sp,  allowMultiple = TRUE)
saveRDS(list(silva = tax_silva, gtdb = tax_gtdb), file.path(out_dir, "asv_newtax.rds"))

## 2. Crosswalk 2022 "Prevotella" -> new genera ----------------------------------------
##  - expectation: most reads of SILVA-138 Prevotella (Prevotella_9 clade) go to
##    Segatella; a second block to Leyella; a small oral fraction stays Prevotella
##  - compare with the provisional crosswalk before trusting any downstream number
cw <- tibble(ASV = colnames(seqtab),
             old_genus   = tax_old$Genus[match(ASV, tax_old$Feature.ID)],
             silva_genus = tax_silva[seq_chr, "Genus"],
             gtdb_genus  = tax_gtdb[seq_chr, "Genus"],
             silva_sp    = tax_silva[seq_chr, "Species"],
             gtdb_sp     = tax_gtdb[seq_chr, "Species"],
             total_reads = colSums(seqtab))
prev_cw <- cw %>% filter(old_genus %in% "Prevotella") %>%
  group_by(old_genus, silva_genus, gtdb_genus) %>%
  summarise(n_asv = n(), reads = sum(total_reads), .groups = "drop") %>% arrange(-reads)
print(prev_cw, n = 50)
write.csv(prev_cw, file.path(out_dir, "prevotella_crosswalk_dada2.csv"), row.names = FALSE)
prov <- read.csv("5-Reanalysis-2026/1-input/prevotella_reassignment.csv")
agree <- cw %>% inner_join(prov %>% select(ASV, provisional = call), by = "ASV") %>%
  mutate(provisional = sub("_lowconf", "", provisional))
print(with(agree, table(provisional, gtdb_genus)))

## 3. Genus count tables (new and old) for 01_diet_prognosis_reanalysis.R ----------------
mk_genus <- function(g) { g[is.na(g)] <- "Unassigned"; t(rowsum(t(seqtab), g)) }
write.csv(mk_genus(tax_gtdb[seq_chr, "Genus"]),  "5-Reanalysis-2026/1-input/genus_newtax_counts.csv")
write.csv(mk_genus(tax_silva[seq_chr, "Genus"]), "5-Reanalysis-2026/1-input/genus_newtax_silva_counts.csv")
write.csv(mk_genus(cw$old_genus),                "5-Reanalysis-2026/1-input/genus_oldtax_counts.csv")

## 4. Enterotypes after re-assignment (Dirichlet multinomial mixtures) --------------------
##  - the cluster driven by Segatella corresponds to the 2022 "Prevotella enterotype"
if (requireNamespace("DirichletMultinomial", quietly = TRUE)) {
  library(DirichletMultinomial)
  g_new <- mk_genus(tax_gtdb[seq_chr, "Genus"]); g_new <- g_new[, colSums(g_new) / sum(g_new) > 1e-4]
  fits  <- lapply(1:5, function(k) dmn(g_new, k, verbose = FALSE))
  best  <- fits[[which.min(sapply(fits, laplace))]]
  et    <- apply(mixture(best), 1, which.max)
  print(table(et)); print(apply(fitted(best), 2, function(x) names(which.max(x))))
  saveRDS(tibble(NO = rownames(g_new), enterotype = et), file.path(out_dir, "enterotype_newtax.rds"))
  ## agreement with the 2022 k-means enterotype
  d22 <- readRDS("3-Output/01_processed.rds")$d
  print(table(old = d22$Cluster_2[match(rownames(g_new), d22$NO.)], new = et))
}
cat("Done 02 (re-taxonomy)\n")
