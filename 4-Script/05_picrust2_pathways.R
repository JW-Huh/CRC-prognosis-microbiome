################################################################################
##  05_picrust2_pathways.R
##  Inferred metagenome (PICRUSt2, MetaCyc pathways) by microbial hazard score
################################################################################
##  Input  : 1-Input/picrust2/path_abun_unstrat.tsv.gz   (405 pathways x 339 samples)
##           1-Input/picrust2/path_cov_unstrat.tsv.gz    (pathway coverage, 144 pathways)
##           1-Input/picrust2/path_abun_contrib.tsv[.gz] (per-ASV contributions, 873 MB;
##             NOT in the repository - set PICRUST_CONTRIB to its location to rebuild)
##           1-Input/picrust2/derived/*.csv              (compact tables derived from it)
##           3-Output/01_processed.rds, 3-Output/03_risk.rds
##  Output : 3-Output/05_pathway_tests.csv, 2-Figure/Fig4*, FigS5*
##  Paper  : Fig. 4 A-G, Figure S5
##  Checks : 144 pathways tested; 15 differential (9 enriched / 6 depleted in M5-low);
##           PFS log-rank p: PWY-6897 0.0087, HISDEG-PWY 0.032
##  PICRUSt2 (v2.4) was run on the rarefied ASV table with default settings:
##     place_seqs.py -> hsp.py (16S, EC) -> metagenome_pipeline.py --strat_out
##     -> pathway_pipeline.py --per_sequence_contrib  (MetaCyc, minpath)
################################################################################

source("4-Script/00_setup.R")
P <- readRDS("3-Output/01_processed.rds"); d <- P$d; tx <- P$tx; NO. <- P$NO.
risk <- readRDS("3-Output/03_risk.rds")$risk
stopifnot(all(risk$NO. == d$NO.))
d$M5_cat <- risk$M5_cat

## 1. Load pathway tables -------------------------------------------------------------
##  - TRFA (taxon relative function abundance) = taxon_rel_function_abun summed over
##    the ASVs of a sample, i.e. pathway abundance per 100 reads of the community
##  - the stratified contribution file is huge; a compact per-sample TRFA table and
##    the family-level contributions of the 15 selected pathways are stored in
##    1-Input/picrust2/derived/ and rebuilt here only when the raw file is available
contrib_file <- Sys.getenv("PICRUST_CONTRIB", "1-Input/picrust2/path_abun_contrib.tsv.gz")
derived_dir  <- "1-Input/picrust2/derived"; dir.create(derived_dir, showWarnings = FALSE)
f_trfa <- file.path(derived_dir, "path_TRFA_by_sample.csv")
f_fam  <- file.path(derived_dir, "path_TRFA_by_family_selected.csv")

sel_path <- c(   # the 15 pathways reported in Fig. 4 (defined up-front so the derived table is self-contained)
  "CALVIN-PWY" = "Calvin-Benson-Bassham cycle",     "PWY-6897" = "thiamine diphosphate salvage II",
  "THISYN-PWY" = "thiamine diphosphate biosynthesis I", "RIBOSYN2-PWY" = "flavin biosynthesis I",
  "PWY-6703" = "preQ0 biosynthesis",                 "PWY0-862" = "(5Z)-dodecenoate biosynthesis I",
  "PWY-5989" = "stearate biosynthesis II",           "PWY-7664" = "oleate biosynthesis IV",
  "PWY-6282" = "palmitoleate biosynthesis I",
  "GLYCOCAT-PWY" = "glycogen degradation I",         "GLYCOLYSIS" = "glycolysis I",
  "PWY-5484" = "glycolysis II",                      "PWY-7199" = "pyrimidine deoxyribonucleoside salvage",
  "HISDEG-PWY" = "L-histidine degradation I",        "PWY-7237" = "myo-/chiro-/scyllo-inositol degradation")

if (file.exists(contrib_file)) {
  cat("Reading", contrib_file, "(takes a minute) ...\n")
  path <- data.table::fread(contrib_file, select = c("sample", "function", "taxon", "taxon_rel_function_abun")) %>%
    setNames(c("ID", "Function", "ASV", "TRFA")) %>% filter(ID %in% NO.)
  path_sum <- path %>% group_by(ID, Function) %>% summarise(TRFA = sum(TRFA)) %>% ungroup()
  path_fam <- path %>% filter(Function %in% names(sel_path)) %>%
    merge(tx[, c("ASV", "Family")], by = "ASV") %>%
    mutate(Family = ifelse(is.na(Family), "unassigned", Family)) %>%
    group_by(ID, Function, Family) %>% summarise(TRFA = sum(TRFA)) %>% ungroup()
  write.csv(path_sum %>% pivot_wider(names_from = Function, values_from = TRFA, values_fill = 0), f_trfa, row.names = FALSE)
  write.csv(path_fam, f_fam, row.names = FALSE)
  rm(path)
} else {
  cat("Stratified file not found - using derived tables in", derived_dir, "\n")
  path_sum <- read.csv(f_trfa, check.names = FALSE) %>% mutate(ID = as.character(ID)) %>%
    pivot_longer(-ID, names_to = "Function", values_to = "TRFA") %>% filter(TRFA > 0)
  path_fam <- read.csv(f_fam) %>% mutate(ID = as.character(ID))
}
path_sum$ID <- as.character(path_sum$ID)

path_unstrat <- read.delim("1-Input/picrust2/path_abun_unstrat.tsv.gz", check.names = FALSE)
path_cov     <- read.delim("1-Input/picrust2/path_cov_unstrat.tsv.gz",  check.names = FALSE)
cat("Pathways: stratified", n_distinct(path_sum$Function), "/ unstratified", nrow(path_unstrat),
    "/ with coverage", nrow(path_cov), "\n")

## 2. Pathways to test ------------------------------------------------------------------
##  - keep pathways with a coverage estimate (144) and present in >= 10 % of samples
pathways <- intersect(path_unstrat$pathway, path_cov$pathway)
prev_ok  <- path_sum %>% count(Function) %>% filter(n >= N_SAMPLE * 0.1) %>% pull(Function)
pathways <- intersect(pathways, prev_ok)
cat("Pathways tested:", length(pathways), "\n")                        # 144

## 3. M5-low ("Good", n = 74) vs M5-moderate/high ("Bad", n = 259) ----------------------
##  - effect size: difference of mean TRFA (Good - Bad); test: Wilcoxon rank-sum, BH-FDR
##  - selection: FDR < 0.001 (-log10 > 3) and |delta TRFA| >= 5
grp <- data.frame(ID = NO., Group = ifelse(d$M5_cat == "Low", "Good", "Bad"))
trfa_wide <- path_sum %>% filter(Function %in% pathways) %>%
  pivot_wider(names_from = Function, values_from = TRFA, values_fill = 0) %>%
  right_join(data.frame(ID = NO.), by = "ID") %>% mutate(across(-ID, ~ replace_na(.x, 0)))
trfa_wide <- trfa_wide[match(NO., trfa_wide$ID), ]
stopifnot(all(trfa_wide$ID == NO.))

func <- lapply(pathways, function(pw) {
  x <- trfa_wide[[pw]]; g <- grp$Group
  data.frame(Function = pw,
             Good = mean(x[g == "Good"]), Bad = mean(x[g == "Bad"]),
             p_two   = wilcox.test(x[g == "Bad"],  x[g == "Good"])$p.value,
             p_good  = wilcox.test(x[g == "Good"], x[g == "Bad"], alternative = "greater")$p.value,
             p_bad   = wilcox.test(x[g == "Bad"],  x[g == "Good"], alternative = "greater")$p.value)
}) %>% bind_rows() %>%
  mutate(Diff = Good - Bad,
         FDR_two  = -log10(p.adjust(p_two,  "fdr")),
         FDR_good = -log10(p.adjust(p_good, "fdr")),
         FDR_bad  = -log10(p.adjust(p_bad,  "fdr")),
         Type = ifelse(Diff > 0, "Good", "Bad"),
         Sig  = ifelse(FDR_two > 3 & abs(Diff) >= 5, paste0("Sig_", Type), "ns"),
         Description = sel_path[Function]) %>%
  arrange(Sig, -abs(Diff))
write.csv(func, "3-Output/05_pathway_tests.csv", row.names = FALSE)
cat("\nDifferential pathways:\n"); print(func %>% filter(Sig != "ns") %>% select(Function, Description, Good, Bad, Diff, FDR_two, Sig))
stopifnot(setequal(func$Function[func$Sig != "ns"], names(sel_path)))

## Fig. 4A: effect size vs significance
p <- func %>% ggplot(aes(Diff, FDR_two)) +
  geom_hline(yintercept = 3, linetype = "dotted", colour = "gray50") +
  geom_vline(xintercept = c(-5, 0, 5), linetype = "dotted", colour = "gray50") +
  geom_point(shape = 21, size = 3.5, aes(fill = Sig)) +
  geom_text_repel(aes(label = ifelse(Sig == "ns", NA, Function))) +
  scale_fill_manual(values = c(ns = "gray98", Sig_Bad = "red3", Sig_Good = "green4")) +
  scale_y_continuous(breaks = c(0, 3, 6, 9)) + theme_classic() + theme(legend.position = "none", aspect.ratio = 1) +
  labs(x = "delta TRFA (M5-low minus M5-moderate/high)", y = "-log10 FDR")
ggsave("2-Figure/Fig4a_pathway_scatter.svg", p, width = 8, height = 8)

## 4. TRFA and family-level contribution of the 15 pathways (Fig. 4B-C) ----------------
##  - enriched pathways: contributions summed within the Good group
##  - depleted pathways: contributions summed within the Bad group
##  - top five contributing families shown, remainder pooled as "Others"
contrib_panel <- function(type, file) {
  sel <- func %>% filter(Sig == paste0("Sig_", type)) %>% arrange(.data[[type]]) %>% mutate(Function = factor(Function, levels = Function))
  heat <- sel %>% select(Function, Good, Bad) %>% pivot_longer(-Function, names_to = "Group", values_to = "TRFA") %>%
    ggplot(aes(Function, Group, fill = TRFA)) + geom_tile(colour = "white", linewidth = 1.2) +
    scale_fill_gradient(low = "white", high = ifelse(type == "Good", "steelblue", "red4")) +
    theme_classic() + coord_flip() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text.x = element_blank(),
          axis.text.y = element_text(size = rel(1.6), colour = "black"))
  fdr <- sel %>% ggplot(aes(Function, y = "1")) +
    geom_text(aes(label = round(if (type == "Good") FDR_good else FDR_bad, 1))) + theme_classic() + coord_flip() +
    theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank())
  fam <- path_fam %>% filter(Function %in% sel$Function, ID %in% grp$ID[grp$Group == type]) %>%
    group_by(Function, Family) %>% summarise(TRFA = sum(TRFA)) %>% ungroup() %>%
    mutate(Family = gsub("f__", "", Family))
  top5 <- fam %>% group_by(Family) %>% summarise(TRFA = sum(TRFA)) %>% arrange(-TRFA) %>% head(5) %>% pull(Family)
  bars <- fam %>% mutate(Family2 = factor(ifelse(Family %in% top5, Family, "Others"), levels = rev(c(top5, "Others"))),
                         Function = factor(Function, levels = levels(sel$Function))) %>%
    ggplot(aes(Function, TRFA, fill = Family2)) + geom_bar(stat = "identity", position = "fill", colour = "transparent") +
    scale_fill_manual(values = col_family, na.value = "gray60") + coord_flip() + theme_classic() +
    theme(legend.position = "right", axis.title = element_blank(), axis.text = element_blank())
  ggsave(file, ggarrange(heat, fdr, bars, nrow = 1, widths = c(2.5, 0.5, 4)), width = 8, height = 4.3)
  fam %>% group_by(Function) %>% mutate(share = TRFA / sum(TRFA)) %>% arrange(Function, -share) %>% slice_head(n = 3)
}
cat("\nTop contributing families, enriched pathways (in M5-low):\n");  print(contrib_panel("Good", "2-Figure/Fig4b_enriched_pathways_contribution.svg"), n = 30)
cat("\nTop contributing families, depleted pathways (in M5-mod/high):\n"); print(contrib_panel("Bad",  "2-Figure/Fig4c_depleted_pathways_contribution.svg"), n = 30)

## 5. Pathway abundance and survival (Fig. 4D-E, Figure S5) -------------------------------
##  - TRFA dichotomised at the median (>= median = High)
##  - thiamine diphosphate salvage II (PWY-6897) p = 0.0087, L-histidine degradation I
##    (HISDEG-PWY) p = 0.032; the others are marginal at best (Figure S5)
##  - pathway IDs contain "-" which survfit() cannot take as a term: use "_" versions as column names
safe <- function(x) gsub("-", "_", x)
surv_df <- d %>% select(NO., PFS, PFS_time, OS, OS_time, Stage, Stage_1, Recur, Recur_type, PCo1, PCo2, CD8_epithelium) %>%
  bind_cols(trfa_wide[, names(sel_path)] %>% setNames(safe(names(sel_path))))
for (pw in names(sel_path)) surv_df[[paste0(safe(pw), "_1")]] <- ifelse(surv_df[[safe(pw)]] >= median(surv_df[[safe(pw)]]), "High", "Low")
pw_surv <- lapply(names(sel_path), function(pw) {
  good <- func$Type[func$Function == pw] == "Good"
  km <- km_plot(as.formula(sprintf("Surv(PFS_time, PFS) ~ %s_1", safe(pw))), surv_df, c("High", "Low"),
                c(ifelse(good, "green4", "red3"), "gray70"),
                file = sprintf("2-Figure/%s_PFS_%s.svg", ifelse(pw %in% c("PWY-6897", "HISDEG-PWY"), "Fig4d-e", "FigS5"), pw))
  data.frame(Function = pw, Description = sel_path[pw], Type = ifelse(good, "Good", "Bad"), p_logrank_PFS = km$p)
}) %>% bind_rows()
cat("\nPFS by pathway (median split):\n"); print(pw_surv %>% mutate(p_logrank_PFS = signif(p_logrank_PFS, 3)), row.names = FALSE)
write.csv(pw_surv, "3-Output/05_pathway_logrank_PFS.csv", row.names = FALSE)

## 6. The two significant pathways by stage and progression (Fig. 4F) ---------------------
two <- surv_df %>% select(Stage, Stage_1, Recur, PWY_6897, HISDEG_PWY) %>%
  pivot_longer(c(PWY_6897, HISDEG_PWY), names_to = "Pathway", values_to = "TRFA") %>%
  mutate(Stage_1 = factor(Stage_1, levels = c("Developing", "Advanced")), Recur = factor(Recur, levels = c("Y", "N")))
p_stage <- two %>% ggplot(aes(Stage, TRFA)) +
  geom_jitter(shape = 21, size = 2.5, aes(fill = Stage), width = 0.2) + geom_boxplot(outlier.alpha = 0, alpha = 0.7) +
  scale_fill_manual(values = col_stage) + facet_wrap(~ Pathway, scales = "free") +
  stat_compare_means(comparisons = list(c("1", "4")), method = "wilcox") + theme_classic() +
  theme(legend.position = "none", axis.title = element_blank(), strip.text = element_text(size = rel(1.5)))
p_recur <- two %>% ggplot(aes(Recur, TRFA, fill = Recur)) + geom_boxplot(outlier.alpha = 0, alpha = 0.7) +
  scale_fill_manual(values = col_YN) + facet_grid(Pathway ~ Stage_1, scales = "free") +
  stat_compare_means(comparisons = list(c("Y", "N")), method = "wilcox", tip.length = 0) + theme_classic() +
  theme(legend.position = "none", axis.title = element_blank(), strip.text = element_text(size = rel(1.3)))
ggsave("2-Figure/Fig4f_two_pathways_by_stage_and_progression.svg", ggarrange(p_stage, p_recur, ncol = 2), width = 12, height = 6)

## 7. Correlation with tumour-infiltrating CD8+ T cells (Fig. 4G) -------------------------
##  - Spearman, epithelial CD8 count available for 267 patients
##  - thiamine salvage rho 0.12 (p 0.067); histidine degradation rho -0.15 (p 0.022)
cd8 <- surv_df %>% filter(!is.na(CD8_epithelium))
cor_cd8 <- lapply(c("PWY-6897", "HISDEG-PWY"), function(pw) {
  ct <- cor.test(cd8[[safe(pw)]], cd8$CD8_epithelium, method = "spearman", exact = FALSE)
  data.frame(Function = pw, n = nrow(cd8), rho = unname(ct$estimate), p = ct$p.value)
}) %>% bind_rows()
cat("\nSpearman with CD8+ T cells (epithelium):\n"); print(cor_cd8, row.names = FALSE)
p <- cd8 %>% select(CD8_epithelium, PWY_6897, HISDEG_PWY) %>%
  pivot_longer(-CD8_epithelium, names_to = "Pathway", values_to = "TRFA") %>%
  ggplot(aes(TRFA, CD8_epithelium)) + geom_point(shape = 21, fill = "gray70") + geom_smooth(method = "lm", formula = y ~ x) +
  facet_wrap(~ Pathway, scales = "free_x") + theme_classic() + labs(y = "CD8+ T cells (epithelium)")
ggsave("2-Figure/Fig4g_pathways_vs_CD8.svg", p, width = 8, height = 4)

saveRDS(list(func = func, surv_df = surv_df, pw_surv = pw_surv, cor_cd8 = cor_cd8), "3-Output/05_picrust2.rds")
cat("\nDone 05\n")
