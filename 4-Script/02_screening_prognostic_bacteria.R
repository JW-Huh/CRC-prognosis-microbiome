################################################################################
##  02_screening_prognostic_bacteria.R
##  Enterotype survival  ->  species screening  ->  Monte-Carlo CV  ->  Alistipes
################################################################################
##  Input  : 3-Output/01_processed.rds
##  Output : 3-Output/02_species_screening.csv  (79 species, whole data + CV counts)
##           3-Output/02_screening.rds          (d with Alistipes columns, results)
##           2-Figure/Fig2*, FigS2*, FigS3*
##  Paper  : Fig. 2 A-J, Figure S2, Figure S3
##  Checks : Prevotella-high PFS p = 0.026, OS p = 0.0056
##           adjusted HR  P. piscolens 2.52, D. invisus 2.30, F. nucleatum 2.04,
##           Bacteroides sp. 1.86 ; four species pass FDR < 0.1 & HR > 1.5
################################################################################

source("4-Script/00_setup.R")
P <- readRDS("3-Output/01_processed.rds")
d <- P$d; r <- P$r; tx <- P$tx; list.taxa <- P$list.taxa; NO. <- P$NO.

## 1. Follow-up summary --------------------------------------------------------
##  - PFS: median 827 d overall (885 d non-progressed, 450 d progressed); 64 events
##  - OS : median 43.6 months for survivors; 38 deaths
cat("PFS events:", sum(d$PFS), "; OS events:", sum(d$OS), "\n")
cat("Median PFS_time (d):", median(d$PFS_time), " non-progressed:", median(d$PFS_time[d$PFS == 0]),
    " progressed:", median(d$PFS_time[d$PFS == 1]), "\n")
cat("Median OS_time survivors (months):", round(median(d$OS_time[d$OS == 0]) / 30.4, 1), "\n")

## 2. Enterotype and enterotypic genera (Fig. 2A-D, Figure S2A-B) -----------------
##  - enterotype 1 (Bacteroides) vs 2 (Prevotella): PFS p = 0.05, OS p = 0.15
##  - Prevotella dichotomised at the median (15.8 %): PFS p = 0.026, OS p = 0.0056
##  - Bacteroides / Faecalibacterium: no association
##  - cut-offs on the relative-abundance scale, for the Methods section:
for (g in c("g__Prevotella", "g__Bacteroides", "g__Faecalibacterium"))
  cat(sprintf("%-20s High >= %.2f %% (Low max %.2f %%)\n", g,
              rel_pct(min(d[[g]][d[[paste0(g, "_1")]] == "High"])),
              rel_pct(max(d[[g]][d[[paste0(g, "_1")]] == "Low"]))))

km_plot(Surv(PFS_time, PFS) ~ Cluster_2, d, c("1", "2"), col_enterotype,
        "2-Figure/Fig2a_PFS_enterotype.svg")
km_plot(Surv(OS_time, OS) ~ Cluster_2, d, c("1", "2"), col_enterotype,
        "2-Figure/Fig2b_OS_enterotype.svg")
km_plot(Surv(PFS_time, PFS) ~ g__Prevotella_1, d, c("High", "Low"), c("green4", "gray70"),
        "2-Figure/Fig2c_PFS_Prevotella.svg")
km_plot(Surv(OS_time, OS) ~ g__Prevotella_1, d, c("High", "Low"), c("green4", "gray70"),
        "2-Figure/Fig2d_OS_Prevotella.svg")
for (g in c("Bacteroides", "Faecalibacterium")) {
  km_plot(as.formula(sprintf("Surv(PFS_time, PFS) ~ g__%s_1", g)), d, c("High", "Low"),
          c("red3", "gray70"), sprintf("2-Figure/FigS2_PFS_%s.svg", g))
  km_plot(as.formula(sprintf("Surv(OS_time, OS) ~ g__%s_1", g)), d, c("High", "Low"),
          c("red3", "gray70"), sprintf("2-Figure/FigS2_OS_%s.svg", g))
}
cat("\nLog-rank p (PFS / OS)\n")
for (v in c("Cluster_2", "g__Prevotella_1", "g__Bacteroides_1", "g__Faecalibacterium_1"))
  cat(sprintf("  %-24s %.4f / %.4f\n", v,
              logrank_p(as.formula(paste0("Surv(PFS_time, PFS) ~ ", v)), d),
              logrank_p(as.formula(paste0("Surv(OS_time, OS) ~ ", v)), d)))

## 3. Species screening on the whole data set ------------------------------------
##  - candidates: 79 species with prevalence >= 10 % and abundance >= 0.01 %
##  - abundance dichotomised at the median (for most species = detected vs not)
##  - (a) univariate log-rank test, BH-adjusted
##  - (b) multivariate Cox adjusted for T (0-2 vs 3/4), N (0 vs 1/2), M (0 vs 1),
##        neoadjuvant treatment, postoperative chemotherapy, lymphatic invasion,
##        CEA (> 5 ng/ml); 6 patients without CEA are dropped by coxph()
##  - the first submission adjusted for Stage / CEA / lymphatic invasion only;
##    the full covariate set was requested in revision and is what the paper reports
species_names <- list.taxa$Species %>%
  filter(`Prev%` >= PREV_MIN, Rel.abund >= ABUND_MIN) %>% pull(Taxa) %>% clean_taxon()
stopifnot(length(species_names) == 79, all(species_names %in% colnames(d)))

cov <- d %>% transmute(
  PFS, PFS_time,
  Stage_T    = ifelse(`T` < 3, "0", "1"),
  Stage_N    = ifelse(N  < 1, "0", "1"),
  Stage_M    = ifelse(M == 0, "0", "1"),
  NeoadjTx   = ifelse(NeoadjTx == 0, "0", "1"),
  Post_chemo = ifelse(Postoperative_chemotherapy == 0, "0", "1"),
  lym_inv    = ifelse(`lymphatic invasion` == 0, "0", "1"),
  CEA        = ifelse(CEA > 5, "1", "0"))
cox_covariates <- "Stage_T + Stage_N + Stage_M + NeoadjTx + Post_chemo + lym_inv + CEA"

sp_abund <- d[, species_names]                      # continuous counts (rows = samples)
sp_bin   <- sp_abund %>% mutate_all(~ ifelse(.x > median(.x), "1", "0"))

res_whole <- lapply(species_names, function(s) {
  dd <- cov; dd$Bug <- sp_bin[[s]]
  cx <- cox_row(as.formula(paste("Surv(PFS_time, PFS) ~ Bug +", cox_covariates)), dd)
  data.frame(Taxa = s,
             pval_logrank = logrank_p(Surv(PFS_time, PFS) ~ Bug, dd),
             HR_cox = cx$HR, lower95_cox = cx$lower95, upper95_cox = cx$upper95, pval_cox = cx$p)
}) %>% bind_rows() %>%
  mutate(FDR_logrank = p.adjust(pval_logrank, method = "fdr")) %>%
  merge(list.taxa$Species %>% mutate(Taxa = clean_taxon(Taxa)) %>% select(Taxa, Rel.abund, `Prev%`), by = "Taxa")

## 4. Monte-Carlo cross-validation (1,000 random 70/30 splits) -------------------
##  - discovery set = 233 patients, validation set = the remaining 100
##  - abundance re-dichotomised at the median within each set
##  - log-rank: species with p < 0.05 in discovery are re-tested in validation;
##    "cv_logrank_both" = number of splits significant in both sets
##  - Cox: same idea with the adjusted model; a hit requires p < 0.05 and HR > 1
##    in discovery, then p < 0.05 in validation
##  - takes a few minutes (2 x 79 x 1,000 model fits)
split_size <- floor(nrow(d) * 0.7)                   # 233
sp_df <- cbind(cov, sp_abund)                        # rownames = sample NO.
rownames(sp_df) <- rownames(sp_abund) <- rownames(d) <- NO.
bin_within <- function(df) mutate_at(df, vars(all_of(species_names)), ~ ifelse(.x > median(.x), "1", "0"))

## 4-1. log-rank CV
cv_logrank <- matrix(NA_real_, nrow = length(species_names), ncol = 1000, dimnames = list(species_names, NULL))
set.seed(SEED_CV)
for (i in 1:1000) {
  disc  <- sample_n(sp_df, split_size, replace = FALSE)
  valid <- sp_df[!rownames(sp_df) %in% rownames(disc), ]
  if (i == 1) stopifnot(nrow(valid) == nrow(d) - split_size, length(intersect(rownames(disc), rownames(valid))) == 0)
  disc <- bin_within(disc); valid <- bin_within(valid)
  for (s in species_names) {
    p_disc <- logrank_p(as.formula(paste0("Surv(PFS_time, PFS) ~ ", s)), disc)
    if (p_disc < 0.05) cv_logrank[s, i] <- logrank_p(as.formula(paste0("Surv(PFS_time, PFS) ~ ", s)), valid)
  }
}

## 4-2. Cox CV
cv_cox_p <- cv_cox_hr <- matrix(NA_real_, nrow = length(species_names), ncol = 1000, dimnames = list(species_names, NULL))
cox_fit <- function(s, df) {
  cx <- summary(coxph(as.formula(paste("Surv(PFS_time, PFS) ~", s, "+", cox_covariates)), data = df))
  c(p = cx$coefficients[1, 5], HR = cx$conf.int[1, 1])
}
set.seed(SEED_CV)
for (i in 1:1000) {
  disc  <- sample_n(sp_df, split_size, replace = FALSE)
  valid <- sp_df[!rownames(sp_df) %in% rownames(disc), ]
  disc <- bin_within(disc); valid <- bin_within(valid)
  for (s in species_names) {
    fd <- cox_fit(s, disc)
    if (fd["p"] < 0.05 & fd["HR"] > 1) {
      fv <- cox_fit(s, valid); cv_cox_p[s, i] <- fv["p"]; cv_cox_hr[s, i] <- fv["HR"]
    }
  }
}

res_cv <- data.frame(Taxa = species_names,
                     cv_logrank_disc = rowSums(!is.na(cv_logrank)),         # significant in discovery
                     cv_logrank_both = rowSums(cv_logrank < 0.05, na.rm = TRUE),  # ... and in validation
                     cv_cox_disc     = rowSums(!is.na(cv_cox_p)),
                     cv_cox_both     = rowSums(cv_cox_p < 0.05, na.rm = TRUE))

## 5. Candidate pathobionts (Fig. 2J) --------------------------------------------
##  - selection rule in the paper: HR > 1.5 and log-rank FDR < 0.1
##  - four species: Bacteroides sp., F. nucleatum, D. invisus, P. piscolens
res_all <- merge(res_whole, res_cv, by = "Taxa") %>%
  mutate(Taxa_label = Taxa %>% gsub("s__", "", .) %>% gsub("_", " ", .),
         Sig = ifelse(FDR_logrank < 0.1 & HR_cox > 1.5, "Sig", "ns")) %>%
  arrange(-cv_cox_both)
write.csv(res_all, "3-Output/02_species_screening.csv", row.names = FALSE)
cat("\nSelected pathobionts:\n"); print(res_all %>% filter(Sig == "Sig") %>%
  select(Taxa_label, Rel.abund, `Prev%`, pval_logrank, FDR_logrank, HR_cox, pval_cox, cv_logrank_both, cv_cox_both))

p <- res_all %>% filter(HR_cox >= 1) %>%
  ggplot(aes(HR_cox, -log10(FDR_logrank))) +
  geom_hline(yintercept = 1,   linetype = "dotted") +
  geom_vline(xintercept = 1.5, linetype = "dotted") +
  geom_point(aes(size = cv_logrank_both, colour = cv_cox_both)) +
  geom_text_repel(aes(label = ifelse(Sig == "Sig", Taxa_label, NA)), fontface = "italic", size = 5) +
  scale_colour_gradient(low = "gray80", high = "red2") +
  theme_classic() +
  theme(aspect.ratio = 1, legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank(), legend.box = "vertical",
        axis.text = element_text(size = rel(1.1), colour = "black"),
        axis.title = element_text(size = rel(1.6), colour = "black")) +
  labs(x = "Hazard ratio", y = "Log-rank test (-log10 FDR)")
ggsave("2-Figure/Fig2j_screening_biomarkers.svg", p, width = 6, height = 7.5)

## 6. Kaplan-Meier curves of the four pathobionts (Fig. 2F-I; OS in Figure S2C) ----
pathobionts <- c(Bacteroides_sp = "s__Bacteroides_sp.",
                 F_nucleatum    = "s__Fusobacterium_nucleatum",
                 D_invisus      = "s__Dialister_invisus",
                 P_piscolens    = "s__Pyramidobacter_piscolens")
for (nm in names(pathobionts)) {
  v <- paste0(pathobionts[nm], "_1")               # High = detected, Low = not detected
  km_plot(as.formula(sprintf("Surv(PFS_time, PFS) ~ `%s`", v)), d, c("Detected", "ND"),
          c("red3", "gray70"), sprintf("2-Figure/Fig2f-i_PFS_%s.svg", nm))
  km_plot(as.formula(sprintf("Surv(OS_time, OS) ~ `%s`", v)), d, c("Detected", "ND"),
          c("red3", "gray70"), sprintf("2-Figure/FigS2c_OS_%s.svg", nm))
}

## 7. "Bacteroides sp." is Alistipes (Figure S3) ------------------------------------
##  - SILVA 138 labels 8 ASVs as s__Bacteroides_sp.; 7 of them sit in g__Alistipes
##  - the 7 Alistipes ASVs cover every patient positive for "Bacteroides sp."
##  - a tree against 14 Alistipes type strains (MEGA, 1-Input/tree/) places them
##    with A. dispar (3 ASVs) and A. senegalensis (4 ASVs); IDs fixed from the tree
bac_asv <- tx %>% filter(Species == "s__Bacteroides_sp.") %>% pull(ASV)
bac_tab <- tx %>% filter(ASV %in% bac_asv) %>%
  mutate(Count = rowSums(r[ASV, ]), Prevalence = rowSums(r[ASV, ] > 0)) %>%
  arrange(Genus, -Prevalence)
print(bac_tab %>% select(ASV, Genus, Species, Count, Prevalence))
bac_pos <- NO.[colSums(r[bac_asv, ]) > 0]           # 119 positive patients
cat("Patients positive for Bacteroides sp.:", length(bac_pos), "\n")

tile <- r[bac_asv, bac_pos] %>% rownames_to_column("ASV") %>%
  pivot_longer(-ASV, names_to = "ID", values_to = "Count") %>%
  merge(bac_tab[, c("ASV", "Genus")], by = "ASV") %>%
  mutate(Det = case_when(Count == 0 ~ "ND", Genus == "g__Alistipes" ~ "Alistipes", TRUE ~ "non-Alistipes"),
         ASV = factor(ASV, levels = bac_tab$ASV))
p <- ggplot(tile, aes(ID, ASV, fill = Det)) + geom_tile() +
  scale_fill_manual(values = c("Alistipes" = "red3", "non-Alistipes" = "orange2", "ND" = "gray95")) +
  scale_y_discrete(limits = rev) + theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title = element_blank(), axis.text.y = element_text(colour = "black"))
ggsave("2-Figure/FigS3a_Bacteroides_sp_ASV_detection.svg", p, width = 10, height = 2.5)

p <- ggplot(bac_tab, aes(x = "", y = Count, fill = ASV)) + geom_bar(stat = "identity", width = 1) +
  coord_polar("y") + theme_void() + theme(legend.position = "none")
ggsave("2-Figure/FigS3b_Bacteroides_sp_ASV_abundance.svg", p, width = 4, height = 4)

if (requireNamespace("ggtree", quietly = TRUE)) {
  suppressPackageStartupMessages(library(ggtree))
  tree <- ape::read.tree("1-Input/tree/Alistipes_tree.nwk")
  lab  <- data.frame(label = tree$tip.label) %>%
    mutate(type = ifelse(label %in% bac_asv, "ASV", "type strain"))
  p <- ggtree(tree) %<+% lab +
    geom_tiplab(aes(colour = type), size = 3.5) +
    scale_colour_manual(values = c("ASV" = "red2", "type strain" = "black")) +
    theme(legend.position = "none") + ggplot2::xlim(0, 0.45)
  ggsave("2-Figure/FigS3c_Alistipes_tree.svg", p, width = 10, height = 5)
}

A.dispar       <- c("690acd4dd0fcb54a27f460104371e964", "535ca2f9c39956853b755d1c1b464606",
                    "ccf990c93977826a69bce80a2079f6e3")
A.senegalensis <- c("7f2afdb722ec5a6a4516d0bac56cef72", "31c40d225eae3bc052fe740ad5a2d893",
                    "c997f158362640120951a6caebf93ade", "65ae54e761425549c369601e23f4eba6")
d <- d %>% mutate(
  s__Alistipes_dispar         = colSums(r[A.dispar, ]),
  s__Alistipes_senegalensis   = colSums(r[A.senegalensis, ]),
  s__Alistipes_dispar_1       = ifelse(s__Alistipes_dispar > 0, "High", "Low"),
  s__Alistipes_senegalensis_1 = ifelse(s__Alistipes_senegalensis > 0, "High", "Low"))
km_plot(Surv(PFS_time, PFS) ~ s__Alistipes_dispar_1, d, c("Detected", "ND"), c("red3", "gray70"),
        "2-Figure/FigS3d_PFS_Alistipes_dispar.svg")
km_plot(Surv(PFS_time, PFS) ~ s__Alistipes_senegalensis_1, d, c("Detected", "ND"), c("red3", "gray70"),
        "2-Figure/FigS3e_PFS_Alistipes_senegalensis.svg")
cat("A. dispar PFS p:",       signif(logrank_p(Surv(PFS_time, PFS) ~ s__Alistipes_dispar_1, d), 3),
    " A. senegalensis PFS p:", signif(logrank_p(Surv(PFS_time, PFS) ~ s__Alistipes_senegalensis_1, d), 3), "\n")

## 8. Save --------------------------------------------------------------------------
saveRDS(list(d = d, res_all = res_all, pathobionts = pathobionts,
             cv_logrank = cv_logrank, cv_cox_p = cv_cox_p, cv_cox_hr = cv_cox_hr),
        "3-Output/02_screening.rds")
cat("\nSaved 3-Output/02_screening.rds\n")
