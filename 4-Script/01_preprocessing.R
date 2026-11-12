################################################################################
##  01_preprocessing.R
##  Raw ASV table  ->  rarefied counts  ->  collapsed taxa  ->  master table `d`
################################################################################
##  Input  : 1-Input/CRC_prognosis_Metadata_333.csv   (333 patients x 106 variables)
##           1-Input/CRC_prognosis_ASV_table.csv      (17,194 ASVs x 339 samples)
##           1-Input/CRC_prognosis_Taxonomy.csv       (SILVA 138 assignment per ASV)
##  Output : 3-Output/01_processed.rds   list(m, r, tx, list.taxa, taxa_summary,
##                                            d, crc.d, eig, NO.)
##           3-Output/01_*_counts_rarefied.csv   collapsed count tables
##  Checks : 16,008 ASVs after rarefaction; enterotype 1/2 = 96/237 patients;
##           median Prevotella 15.8 %; PCo1 7.84 %, PCo2 4.28 %
################################################################################

source("4-Script/00_setup.R")

## 1. Import -------------------------------------------------------------------
##  - 339 patients were recruited; 6 were excluded (2 carcinoma in situ,
##    2 refused surgery, 2 non-adenocarcinoma) -> Metadata_333 holds the final set
##  - the ASV table still has 339 columns; keep only the 333 included patients
m   <- read.csv("1-Input/CRC_prognosis_Metadata_333.csv", check.names = FALSE)
NO. <- as.character(m$NO.)                       # sample IDs, used everywhere as keys

a <- read.csv("1-Input/CRC_prognosis_ASV_table.csv", row.names = 1, check.names = FALSE)
a <- a[, NO.]
stopifnot(nrow(a) == 17194, ncol(a) == N_SAMPLE)
cat("Total reads in 333 samples:", sum(a), "\n")   # 34,184,627 reads

## Taxonomy string "d__; p__; c__; o__; f__; g__; s__" -> one column per rank
tx_raw <- read.csv("1-Input/CRC_prognosis_Taxonomy.csv")
ranks  <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
tx <- strsplit(tx_raw$Taxon, "; ") %>%
  lapply(function(x) c(x, rep(NA, 7 - length(x)))) %>%
  do.call(rbind, .) %>% as.data.frame() %>% set_colnames(ranks) %>%
  mutate(ASV = tx_raw$Feature.ID, Confidence = tx_raw$Confidence) %>%
  select(ASV, Phylum, Family, Genus, Species)

## 2. Sequencing depth ---------------------------------------------------------
depth <- colSums(a)
cat("Depth: min", min(depth), "/ median", median(depth), "/ max", max(depth), "\n")
##  - min 64,604 ; median 102,049 ; max 153,429  -> rarefy to the minimum

## 3. Rarefaction ----------------------------------------------------------------
##  - even library size: every sample is subsampled to 64,604 reads without replacement
##  - seed fixed; vegan 2.5-7 (2021) and 2.7-6 (2026) give the identical table
##  - ASVs that lose every read are dropped: 17,194 -> 16,008 ASVs
set.seed(SEED_RAREFY)
r <- a %>% t() %>% rrarefy(sample = RAREFY_DEPTH) %>% t() %>%
  as.data.frame() %>% set_colnames(NO.) %>% filter(rowSums(.) > 0)
cat("ASVs after rarefaction:", nrow(r), "\n")
rm(a)

## 4. Collapse ASVs to taxonomic ranks ---------------------------------------------
##  - counts   : taxon x sample table (NA rank -> "unassigned")
##  - summary  : total count, relative abundance (% of all reads), prevalence
##  - the summary drops uninformative labels (unassigned / uncultured / metagenome ...)
collapse_rank <- function(r, tx, rank, drop = c("unassigned", "uncultured")) {
  counts <- r %>% rownames_to_column("ASV") %>%
    merge(tx[, c("ASV", rank)], by = "ASV") %>%
    select(all_of(c(rank, NO.))) %>%
    group_by(.data[[rank]]) %>% summarise_all(sum) %>% as.data.frame()
  counts[[rank]] <- ifelse(is.na(counts[[rank]]), "unassigned", counts[[rank]])
  counts <- column_to_rownames(counts, rank)
  summ <- data.frame(Count = rowSums(counts)) %>%
    mutate(Rel.abund  = round(Count / (N_SAMPLE * RAREFY_DEPTH) * 100, 5),
           Prevalence = rowSums(counts > 0),
           `Prev%`    = round(Prevalence / N_SAMPLE * 100, 1)) %>%
    arrange(-Rel.abund, -`Prev%`) %>% rownames_to_column("Taxa") %>%
    filter(!grepl(paste(drop, collapse = "|"), Taxa))
  list(counts = counts, summary = summ)
}
##  - species labels that carry no taxonomic information are dropped as well
##    ("uncultured_bacterium", "gut_metagenome", "human_gut")
##  - NOTE "s__unidentified" (a SILVA placeholder) is deliberately kept: it slipped
##    through the published filter (typo in the pattern) and is one of the 79
##    screened species, so keeping it preserves the FDR values reported in the paper
phylum  <- collapse_rank(r, tx, "Phylum")    # Bacteroidota 47.6 %, Firmicutes 42.1 %
family  <- collapse_rank(r, tx, "Family")    # Prevotellaceae 21.3 %, Bacteroidaceae 19.6 %
genus   <- collapse_rank(r, tx, "Genus")     # Bacteroides 19.6 %, Prevotella 19.4 %
species <- collapse_rank(r, tx, "Species",   # B. vulgatus 5.4 %, B. plebeius 4.4 %, P. copri 2.3 %
                         drop = c("unassigned", "uncultured", "metagenome", "_bacterium", "_gut"))

list.taxa <- list(p = phylum$counts,  Phylum  = phylum$summary,
                  f = family$counts,  Family  = family$summary,
                  g = genus$counts,   Genus   = genus$summary,
                  s = species$counts, Species = species$summary)
for (k in c("Phylum", "Family", "Genus", "Species")) { cat("\nTop 5", k, "\n"); print(head(list.taxa[[k]], 5)) }

## Long summary table used for ranking / barplots in later scripts
taxa_summary <- bind_rows(
  list.taxa$Family  %>% mutate(Type = "Family",  Rank = row_number()),
  list.taxa$Genus   %>% mutate(Type = "Genus",   Rank = row_number()),
  list.taxa$Species %>% mutate(Type = "Species", Rank = row_number()))

## Export collapsed tables (handy for re-use outside R)
write.csv(genus$counts,   "3-Output/01_genus_counts_rarefied.csv")
write.csv(species$counts, "3-Output/01_species_counts_rarefied.csv")
write.csv(family$counts,  "3-Output/01_family_counts_rarefied.csv")

## 5. Taxa to be tested for prognosis -------------------------------------------
##  - genus + species with prevalence >= 10 % and relative abundance >= 0.01 %
##  - 209 taxa in total (130 genera + 79 species); only the 79 species enter the
##    final biomarker screen (script 02), genera are kept for exploration
##  - "_1" suffix = dichotomised at the median (High / Low)
taxa <- rbind(
  genus$counts[  list.taxa$Genus   %>% filter(`Prev%` >= PREV_MIN, Rel.abund >= ABUND_MIN) %>% pull(Taxa), ],
  species$counts[list.taxa$Species %>% filter(`Prev%` >= PREV_MIN, Rel.abund >= ABUND_MIN) %>% pull(Taxa), ]) %>%
  t() %>% as.data.frame()
colnames(taxa) <- clean_taxon(colnames(taxa))
taxa_1 <- taxa %>% apply(2, hi_lo) %>% as.data.frame() %>% set_colnames(paste0(colnames(taxa), "_1"))
taxa_to_test   <- colnames(taxa)
taxa_to_test_1 <- colnames(taxa_1)
cat("\nTaxa to test:", length(taxa_to_test), "\n")

## 6. Master table `d` = metadata + taxa + derived host variables -----------------
##  - blood inflammation indices with published cut-offs (see Methods)
##  - Stage_1: 0-2 "Developing" vs 3-4 "Advanced"
##  - CEA_1 in the metadata is coded 1 (<= 5 ng/ml) / 2 (> 5 ng/ml)
d <- m %>%
  merge(cbind(taxa, taxa_1) %>% rownames_to_column("NO."), by = "NO.") %>%
  mutate(
    PB_ratio  = log2((g__Prevotella + 1) / (g__Bacteroides + 1)),
    NLR = neutrophil / lymphocyte,          NLR_1 = ifelse(NLR >= 5,   "High", "Low"),
    PLR = platelet   / lymphocyte,          PLR_1 = ifelse(PLR >= 150, "High", "Low"),
    LMR = lymphocyte / monocyte,            LMR_1 = ifelse(LMR >= 2.4, "High", "Low"),
    PNI = lymphocyte * 0.005 + Albumin * 10, PNI_1 = ifelse(PNI >= 45,  "High", "Low"),
    Stage     = as.character(Stage),
    Stage_num = as.numeric(Stage),
    Stage_1   = ifelse(Stage_num >= 3, "Advanced", "Developing"),
    Shannon      = diversity(t(r)),
    Observed_ASV = colSums(r > 0),
    CEA_1        = ifelse(CEA_1 > 1, "High", "Low"),
    fibrinogen_1 = ifelse(fibrinogen >= 335, "High", "Low"),
    Age_1        = ifelse(Age >= 65, "High", "Low"),
    BMI_1        = ifelse(BMI >= 25, "High", "Low"),
    preoperative_medication = as.character(preoperative_medication),
    PFS_time = as.numeric(PFS_time))
rm(taxa, taxa_1)
stopifnot(all(d$NO. == NO.))                      # merge() keeps the original order here
cat("Median CEA", median(m$CEA, na.rm = TRUE), "; CEA-high n =", sum(d$CEA_1 == "High", na.rm = TRUE), "\n")

## 7. Beta diversity: Bray-Curtis -> PCoA --------------------------------------
crc.d <- vegdist(t(r))                            # Bray-Curtis on rarefied ASV counts
pcoa  <- cmdscale(crc.d, k = 10, eig = TRUE)
d     <- bind_cols(as.data.frame(pcoa$points) %>% set_colnames(paste0("PCo", 1:10)), d)
eig   <- pcoa$eig; eig <- round(eig[eig > 0] / sum(eig[eig > 0]) * 100, 2)[1:10]
cat("Variance explained: PCo1", eig[1], "% ; PCo2", eig[2], "%\n")   # 7.84 / 4.28

## 8. Enterotype by k-means (k = 2) ----------------------------------------------
##  - k-means applied to the Bray-Curtis distance matrix (each sample = its row of
##    distances); k = 2 chosen heuristically; nstart = 100 to avoid local optima
##  - labels made deterministic: "1" = Bacteroides-type, "2" = Prevotella-type
set.seed(SEED_KMEANS)
km <- kmeans(crc.d, centers = 2, iter.max = 10000, nstart = 100)$cluster
prev_by_cl <- tapply(d$g__Prevotella, km, mean)
d$Cluster_2 <- ifelse(km == names(which.max(prev_by_cl)), "2", "1")
print(table(d$Cluster_2))                         # 1 = 96, 2 = 237

## 9. Progression categories -------------------------------------------------------
##  - Recur     : any progression during follow-up (PFS event, 64 / 333)
##  - Recur_1y/2y/3y : early recurrence within 1/2/3 years; censored before the
##                     landmark -> "NA" (not evaluable)
##  - Recur_type: "yet" (no event), "<1y" (within 1 year), ">1y" (later)
yr <- 365
d <- d %>% mutate(
  Recur    = ifelse(PFS == 1, "Y", "N"),
  Recur_1y = case_when(PFS == 1 & PFS_time <= yr     ~ "Y", PFS == 1 ~ "N",
                       PFS == 0 & PFS_time <= yr     ~ "NA", TRUE ~ "N"),
  Recur_2y = case_when(PFS == 1 & PFS_time <= 2 * yr ~ "Y", PFS == 1 ~ "N",
                       PFS == 0 & PFS_time <= 2 * yr ~ "NA", TRUE ~ "N"),
  Recur_3y = case_when(PFS == 1 & PFS_time <= 3 * yr ~ "Y", PFS == 1 ~ "N",
                       PFS == 0 & PFS_time <= 3 * yr ~ "NA", TRUE ~ "N"),
  Recur_type = case_when(PFS == 0 ~ "yet", Recur_1y == "Y" ~ "<1y", TRUE ~ ">1y"),
  Recur_type = factor(Recur_type, levels = c("yet", ">1y", "<1y")))
print(with(d, table(Recur)))                      # 64 progressed

## 10. Save ------------------------------------------------------------------------
saveRDS(list(m = m, r = r, tx = tx, list.taxa = list.taxa, taxa_summary = taxa_summary,
             taxa_to_test = taxa_to_test, d = d, crc.d = crc.d, eig = eig, NO. = NO.),
        "3-Output/01_processed.rds")
cat("\nSaved 3-Output/01_processed.rds\n")
