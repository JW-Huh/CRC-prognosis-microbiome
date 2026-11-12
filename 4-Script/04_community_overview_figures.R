################################################################################
##  04_community_overview_figures.R
##  Community-level description: taxonomy barplots, PCoA / PERMANOVA, alpha
##  diversity, enterotype vs. clinical variables, host prognostic factors
################################################################################
##  Input  : 3-Output/01_processed.rds, 3-Output/03_risk.rds
##  Output : 3-Output/04_permanova.csv, 04_enterotype_vs_clinical.csv,
##           04_host_variables_logrank.csv
##           2-Figure/Fig1*, FigS1*, FigS4*
##  Paper  : Fig. 1 B-H, Figure S1, Figure S4
##  Checks : PERMANOVA stage p = 0.012, progression p = 0.011, death p = 0.007
##           T p = 0.012-0.016, N p = 0.063, M p = 0.01 (999 permutations, seed 1)
##           Wilcoxon by progression: Akkermansiaceae p = 0.0039, Prevotellaceae 0.0228
################################################################################

source("4-Script/00_setup.R")
P <- readRDS("3-Output/01_processed.rds")
d <- P$d; r <- P$r; list.taxa <- P$list.taxa; taxa_summary <- P$taxa_summary; crc.d <- P$crc.d; eig <- P$eig
NO. <- P$NO.
d$M5_cat <- readRDS("3-Output/03_risk.rds")$risk$M5_cat

## 1. Family-level taxonomic profile (Fig. 1B, 1G) ---------------------------------
##  - six most abundant families + Fusobacteriaceae + Akkermansiaceae, rest = Other
##  - samples ordered by Bacteroidaceae (descending) then Prevotellaceae
fam_show <- c(taxa_summary %>% filter(Type == "Family", Rank <= 6) %>% pull(Taxa),
              "f__Fusobacteriaceae", "f__Akkermansiaceae")
fam_rel <- list.taxa$f[fam_show, NO.] %>% t() %>% as.data.frame() %>% `/`(RAREFY_DEPTH) %>%
  mutate(Other = 1 - rowSums(.)) %>% rownames_to_column("NO.") %>%
  merge(d %>% select(NO., Cluster_2, Stage, Stage_1, Recur, OS), by = "NO.")
fam_levels <- c("Other", gsub("f__", "", rev(fam_show)))
fam_long <- function(df, group) df %>%
  group_by(.data[[group]]) %>% summarise(across(all_of(c(fam_show, "Other")), mean)) %>%
  pivot_longer(-all_of(group), names_to = "Family", values_to = "Rel.abund") %>%
  mutate(Family = factor(gsub("f__", "", Family), levels = fam_levels))

sample_order <- fam_rel %>% arrange(-f__Bacteroidaceae, f__Prevotellaceae) %>% pull(NO.)
p <- fam_rel %>% select(NO., all_of(fam_show), Other) %>%
  pivot_longer(-NO., names_to = "Family", values_to = "Rel.abund") %>%
  mutate(Family = factor(gsub("f__", "", Family), levels = fam_levels),
         NO. = factor(NO., levels = sample_order)) %>%
  ggplot(aes(NO., Rel.abund, fill = Family)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = col_family) + theme_classic() +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks.x = element_blank(), aspect.ratio = 0.8)
ggsave("2-Figure/Fig1b_family_barplot_per_patient.svg", p, width = 6, height = 5.5)

## mean composition by enterotype / stage / progression / death (Fig. 1G)
fam_rel <- fam_rel %>% mutate(Stage_4 = factor(ifelse(Stage %in% c("0", "1"), "0/1", Stage), levels = c("0/1", "2", "3", "4")),
                              Death = ifelse(OS == 1, "Y", "N"))
bar_by <- function(group, ratio) fam_long(fam_rel, group) %>%
  ggplot(aes(.data[[group]], Rel.abund, fill = Family)) + geom_bar(stat = "identity") +
  scale_fill_manual(values = col_family) + theme_classic() +
  theme(legend.position = "none", axis.title = element_blank(), axis.text = element_blank(), aspect.ratio = ratio)
ggsave("2-Figure/Fig1g_family_by_enterotype_stage_progression_death.svg",
       ggarrange(bar_by("Cluster_2", 6), bar_by("Stage_4", 3), bar_by("Recur", 6), bar_by("Death", 6),
                 nrow = 1, widths = c(1, 1.79, 1, 1)), width = 8, height = 5.5)

## Wilcoxon tests quoted in the text (family abundance by progression / enterotype)
fam_test <- expand_grid(Family = c("f__Akkermansiaceae", "f__Fusobacteriaceae", "f__Prevotellaceae"),
                        Group = c("Recur", "Cluster_2")) %>%
  rowwise() %>%
  mutate(p = wilcox.test(fam_rel[[Family]] ~ fam_rel[[Group]])$p.value) %>% ungroup()
cat("\nFamily abundance, Wilcoxon:\n"); print(fam_test %>% mutate(p = signif(p, 3)))

## 2. Enterotype on the PCoA plane (Fig. 1C) ------------------------------------------
##  - 75 % data ellipse per enterotype, spokes to the robust centroid (MASS::cov.trob)
##  - marginal loess curves of Prevotella / Bacteroides counts along PCo1 and PCo2
cen <- d %>% group_by(Cluster_2) %>%
  summarise(cx = MASS::cov.trob(cbind(PCo1, PCo2))$center[1], cy = MASS::cov.trob(cbind(PCo1, PCo2))$center[2])
p_main <- d %>% left_join(cen, by = "Cluster_2") %>%
  ggplot(aes(PCo1, PCo2, colour = Cluster_2, fill = Cluster_2)) +
  geom_segment(aes(xend = cx, yend = cy), alpha = 0.3, linewidth = 0.3) +
  geom_point(shape = 21, size = 2.5, alpha = 0.8) +
  stat_ellipse(type = "t", level = 0.75, geom = "polygon", alpha = 0.2, colour = NA) +
  geom_label(data = cen, aes(cx, cy, label = Cluster_2), colour = "black", fill = "white") +
  scale_fill_manual(values = alpha(c("pink", "lightgreen"), 1)) +
  scale_colour_manual(values = col_enterotype) + theme_classic() +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title.x = element_text(size = rel(1.6), vjust = -1), axis.title.y = element_text(size = rel(1.6), vjust = 2)) +
  labs(x = paste0("PCo1 (", eig[1], "%)"), y = paste0("PCo2 (", eig[2], "%)"))
marg <- function(axis) d %>% select(all_of(axis), g__Prevotella, g__Bacteroides) %>%
  pivot_longer(-all_of(axis), names_to = "Genus", values_to = "Count") %>%
  ggplot(aes(.data[[axis]], Count, colour = Genus, fill = Genus, group = Genus)) +
  geom_smooth(alpha = 0.3, method = "loess", span = 0.9, formula = y ~ x) +
  scale_fill_manual(values = c(g__Bacteroides = "pink", g__Prevotella = "lightgreen")) +
  scale_colour_manual(values = c(g__Bacteroides = "red3", g__Prevotella = "green4")) +
  theme_classic() + theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
ggsave("2-Figure/Fig1c_PCoA_enterotype.svg",
       ggarrange(marg("PCo1"), NULL, p_main, marg("PCo2") + rotate(), ncol = 2, nrow = 2,
                 align = "hv", widths = c(3, 1), heights = c(1, 3)), width = 8, height = 8)

## 3. PERMANOVA on Bray-Curtis (Fig. 1D-F, Figure S1 and text) -----------------------
##  - every clinical variable in the list is tested with the same seed (999 perms)
##  - progression / death / stage: p = 0.011 / 0.007 / 0.012
d <- d %>% mutate(PFS_c = as.character(PFS), OS_c = as.character(OS), T_c = as.character(`T`),
                  N_c = as.character(N), M_c = as.character(M),
                  Stage_01 = ifelse(Stage == "0", "1", Stage),           # stage 0 pooled with 1 (7 patients)
                  NeoadjTx_1 = ifelse(NeoadjTx == 0, "0", "1"),
                  Post_chemo = as.character(Postoperative_chemotherapy),
                  Comorbidity = as.character(Comorbidity),
                  Lym_inv = as.character(`lymphatic invasion`),
                  Response_chemo2 = ifelse(Postoperative_chemotherapy == 1,
                                           ifelse(Response_chemo %in% c("CR", "PR", "SD", "Mixed response"), "responder", "non"), NA),
                  across(c(`antihypertensive drug`, `lipid-lowering drug`, `antidiabetic drug`,
                           `Nonsteroidal anti-inflammatory drug`, `antiplatelet drug`, `antiulcer drug`,
                           `drug for benign prostatic hypertrophy`), ~ ifelse(is.na(.x), "0", "1"), .names = "{.col}_01"),
                  across(c(Hypertension, Diabetes, `Heart disease`, `Pulmonary disease`, `Liver disease`,
                           `Cerebrovascular accident`, `Renal disease`, `Thyroid disease`), as.character, .names = "{.col}_c"),
                  Recur_1y_v = na_if(Recur_1y, "NA"), Recur_2y_v = na_if(Recur_2y, "NA"))
perm_vars <- c("Cluster_2", "Stage", "Stage_01", "Stage_1", "PFS_c", "OS_c", "Recur_type", "Recur_1y_v", "Recur_2y_v",
               "T_c", "N_c", "M_c", "Sex", "CEA_1", "fibrinogen_1", "Lym_inv", "NeoadjTx_1", "Post_chemo",
               "Response_chemo2", "Comorbidity", "preoperative_medication",
               paste0(c("antihypertensive drug", "lipid-lowering drug", "antidiabetic drug",
                        "Nonsteroidal anti-inflammatory drug", "antiplatelet drug", "antiulcer drug",
                        "drug for benign prostatic hypertrophy"), "_01"),
               paste0(c("Hypertension", "Diabetes", "Heart disease", "Pulmonary disease", "Liver disease",
                        "Cerebrovascular accident", "Renal disease", "Thyroid disease"), "_c"))
permanova <- data.frame(Variable = perm_vars,
                        n = sapply(perm_vars, function(v) sum(!is.na(d[[v]]))),
                        p = sapply(perm_vars, function(v) permanova_p(crc.d, d, v)))
write.csv(permanova, "3-Output/04_permanova.csv", row.names = FALSE)
cat("\nPERMANOVA (Bray-Curtis, 999 permutations):\n"); print(permanova, row.names = FALSE)
pp <- function(v) permanova$p[permanova$Variable == v]

## PCoA panels coloured by stage / progression / death with marginal densities (Fig. 1D-F)
pcoa_panel <- function(var, cols, p, file) {
  d1 <- d %>% filter(!is.na(.data[[var]]))
  main <- draw_pcoa(d1, .data[[var]], p) + scale_fill_manual(values = cols) + t_none
  dens <- function(axis) ggplot(d1, aes(.data[[axis]], colour = .data[[var]], fill = .data[[var]])) +
    geom_density(alpha = 0.3) + scale_fill_manual(values = cols) + scale_colour_manual(values = cols) +
    theme_void() + theme(legend.position = "none")
  ggsave(file, ggarrange(dens("PCo1"), NULL, main, dens("PCo2") + rotate(), ncol = 2, nrow = 2,
                         align = "hv", widths = c(3, 1), heights = c(1, 3)), width = 8, height = 8)
}
pcoa_panel("Stage_01", col_stage[-1], pp("Stage"), "2-Figure/Fig1d_PCoA_stage.svg")
pcoa_panel("PFS_c", col_01, pp("PFS_c"), "2-Figure/Fig1e_PCoA_progression.svg")
pcoa_panel("OS_c", col_01, pp("OS_c"), "2-Figure/Fig1f_PCoA_death.svg")
pcoa_panel("Recur_type", col_recur, pp("Recur_type"), "2-Figure/FigS1_PCoA_recurrence_type.svg")
ggsave("2-Figure/FigS1_PCoA_CEA.svg", draw_pcoa(d, CEA_1, pp("CEA_1")) + scale_fill_manual(values = c(High = "coral3", Low = "cornsilk2")), width = 5, height = 5.2)
ggsave("2-Figure/FigS1_PCoA_fibrinogen.svg", draw_pcoa(d, fibrinogen_1, pp("fibrinogen_1")) + scale_fill_manual(values = c(High = "coral3", Low = "cornsilk2")), width = 5, height = 5.2)

## 4. Alpha diversity (Fig. 1H) ----------------------------------------------------------
##  - Shannon index on rarefied ASV counts; higher in enterotype 2, otherwise flat
shan_box <- function(var, cols, comparisons, label_y) d %>% filter(!is.na(.data[[var]])) %>%
  ggplot(aes(.data[[var]], Shannon)) +
  geom_violin(fill = "gray90", trim = FALSE, alpha = 0.5, colour = NA) +
  geom_boxplot(outlier.alpha = 0, width = 0.3, aes(fill = .data[[var]])) +
  stat_compare_means(comparisons = comparisons, tip.length = 0, label.y = label_y, method = "wilcox") +
  scale_fill_manual(values = cols) + theme_classic() +
  theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
d <- d %>% mutate(Stage_4 = factor(ifelse(Stage %in% c("0", "1"), "0/1", Stage), levels = c("0/1", "2", "3", "4")),
                  Death = ifelse(OS == 1, "Y", "N"))
ggsave("2-Figure/Fig1h_Shannon_by_enterotype_stage_progression_death.svg",
       ggarrange(shan_box("Cluster_2", col_enterotype, list(c("1", "2")), 5.5),
                 shan_box("Stage_4", c("0/1" = "cornsilk1", "2" = "cornsilk3", "3" = "coral1", "4" = "coral4"),
                          list(c("2", "3"), c("2", "4"), c("0/1", "4")), c(5.3, 5.5, 5.7)),
                 shan_box("Recur", col_YN, list(c("Y", "N")), 5.5),
                 shan_box("Death", col_YN, list(c("Y", "N")), 5.5),
                 nrow = 1, widths = c(1, 1.79, 1, 1)), width = 8, height = 5)

## 5. Enterotype vs. clinical variables (Figure S1 A-I) -----------------------------
##  - categorical: Fisher's exact test + waffle (one tile per patient)
##  - blood indicators: Wilcoxon rank-sum test (CEA and fibrinogen are higher in type 1)
et_tab <- function(var) {
  tab <- table(d$Cluster_2, d[[var]])
  data.frame(Variable = var, p_fisher = fisher.test(tab)$p.value, p_chisq = chisq.test(tab)$p.value)
}
et_cat <- bind_rows(et_tab("Stage_1"), et_tab("PFS"), et_tab("OS"), et_tab("Sex"))
write.csv(et_cat, "3-Output/04_enterotype_vs_clinical.csv", row.names = FALSE)
cat("\nEnterotype vs clinical (Fisher / chi-square):\n"); print(et_cat, row.names = FALSE)
waffle_pair <- function(var, levels, cols, file) {
  w1 <- waffle_plot(table(factor(d[[var]][d$Cluster_2 == "1"], levels = levels)), cols, rows = 5)
  w2 <- waffle_plot(table(factor(d[[var]][d$Cluster_2 == "2"], levels = levels)), cols, rows = 12)
  ggsave(file, ggarrange(w1, w2, nrow = 1, widths = c(5, 12)), width = 6, height = 6)
}
waffle_pair("Stage", as.character(0:4), unname(col_stage), "2-Figure/FigS1a_waffle_stage_by_enterotype.svg")
waffle_pair("PFS", c("0", "1"), c("cornsilk2", "coral3"), "2-Figure/FigS1b_waffle_progression_by_enterotype.svg")
waffle_pair("OS", c("0", "1"), c("cornsilk2", "coral3"), "2-Figure/FigS1c_waffle_death_by_enterotype.svg")
waffle_pair("Sex", c("F", "M"), c("orange2", "steelblue"), "2-Figure/FigS1d_waffle_sex_by_enterotype.svg")

ind_all <- c("lymphocyte", "WBC", "neutrophil", "eosinophil", "basophil", "monocyte",
             "RBC", "HB", "Hct", "MCV", "MCH", "MCHC",
             "platelet", "fibrinogen", "PT_INR", "PT %", "PT sec", "aPTT",
             "glucose", "Cholesterol", "Albumin", "Protein_total", "Bilirubin_total", "Alkaline phosphatase",
             "creatinine", "CRP", "GOT_AST", "Calcium", "Phosphorus", "uric acid",
             "CEA", "NLR", "LMR", "PLR", "PNI", "Shannon")
ind_log <- c("CRP", "CEA")                                  # right-skewed -> log2(x + 1) for display
ind_long <- d %>% select(Cluster_2, all_of(ind_all)) %>%
  pivot_longer(-Cluster_2, names_to = "variable", values_to = "value") %>% filter(!is.na(value)) %>%
  mutate(value = ifelse(variable %in% ind_log, log2(value + 1), value), variable = factor(variable, levels = ind_all))
ind_test <- ind_long %>% group_by(variable) %>% summarise(p_wilcox = wilcox.test(value ~ Cluster_2)$p.value) %>%
  arrange(p_wilcox)
write.csv(ind_test, "3-Output/04_blood_indicators_by_enterotype.csv", row.names = FALSE)
cat("\nBlood indicators by enterotype (lowest p):\n"); print(head(ind_test, 8))
p <- ind_long %>% ggplot(aes(Cluster_2, value)) +
  geom_violin(fill = "gray90", colour = "transparent", alpha = 0.8, trim = FALSE) +
  geom_boxplot(aes(fill = Cluster_2), width = 0.3, outlier.alpha = 0) +
  scale_fill_manual(values = col_enterotype) + facet_wrap(~ variable, scales = "free", nrow = 4) +
  stat_compare_means(comparisons = list(c("1", "2")), size = 3, tip.length = 0, method = "wilcox") +
  theme_classic() + theme(legend.position = "top", axis.title = element_blank(),
                          axis.text.x = element_text(size = rel(0.9), colour = "black"), axis.text.y = element_text(size = rel(0.7))) +
  labs(fill = "Enterotype")
ggsave("2-Figure/FigS1e-i_indicators_by_enterotype.svg", p, width = 14, height = 12)

## 6. Host prognostic factors (Figure S4) -------------------------------------------------
##  - log-rank test for PFS; published cut-offs where they exist, otherwise the median
##  - confirmed: age, stage, lymphatic invasion, CEA, fibrinogen, bilirubin, CRP
##  - not associated: NLR, LMR, PNI, PLR
host_spec <- tribble(
  ~var,                 ~group,        ~labels,                    ~cols,
  "Age_1",              NA,            c(">=65", "<65"),           c("red3", "gray70"),
  "Sex",                NA,            c("F", "M"),                c("orange2", "steelblue"),
  "BMI_1",              NA,            c(">=25", "<25"),           c("black", "gray70"),
  "Stage_1",            NA,            c("3/4", "0/1/2"),          c("red3", "gray70"),
  "CEA_1",              NA,            c("> 5", "<= 5"),           c("red3", "gray70"),
  "Lym_inv",            NA,            c("No", "Yes"),             c("gray70", "red3"),
  "fibrinogen_1",       NA,            c(">=335", "<335"),         c("red3", "gray70"),
  "CRP",                "median",      c("Low", "High"),           c("gray70", "red3"),
  "Alkaline phosphatase", "median",    c("Low", "High"),           c("gray70", "red3"),
  "Albumin",            "median",      c("Low", "High"),           c("gray70", "green4"),
  "Bilirubin_total",    "median",      c("Low", "High"),           c("gray70", "green4"),
  "BUN",                "median",      c("Low", "High"),           c("gray70", "green4"),
  "lymphocyte",         "median",      c("Low", "High"),           c("gray70", "black"),
  "glucose",            "median",      c("Low", "High"),           c("gray70", "black"),
  "Cholesterol",        "median",      c("Low", "High"),           c("gray70", "black"),
  "uric acid",          "median",      c("Low", "High"),           c("gray70", "black"),
  "Protein_total",      "median",      c("Low", "High"),           c("gray70", "black"),
  "creatinine",         "median",      c("Low", "High"),           c("gray70", "black"),
  "Calcium",            "median",      c("Low", "High"),           c("gray70", "black"),
  "Phosphorus",         "median",      c("Low", "High"),           c("gray70", "black"),
  "PNI_1",              NA,            c(">=45", "<45"),           c("green4", "gray70"),
  "NLR_1",              NA,            c(">=5", "<5"),             c("black", "gray70"),
  "LMR_1",              NA,            c(">=2.4", "<2.4"),         c("black", "gray70"),
  "PLR_1",              NA,            c(">=150", "<150"),         c("black", "gray70"))
host_res <- list(); host_panels <- list()
for (i in seq_len(nrow(host_spec))) {
  v <- host_spec$var[i]
  dd <- d %>% filter(!is.na(.data[[v]]))
  ## median-split variables: factor Low < High so that labels and colours line up;
  ## pre-coded variables keep their alphabetical strata order (High/Low, 0/1, F/M ...)
  dd$grp <- if (identical(host_spec$group[[i]], "median"))
    factor(ifelse(dd[[v]] >= median(dd[[v]]), "High", "Low"), levels = c("Low", "High")) else as.character(dd[[v]])
  km <- km_plot(Surv(PFS_time, PFS) ~ grp, dd, host_spec$labels[[i]], host_spec$cols[[i]])
  host_res[[v]]    <- data.frame(Variable = v, n = nrow(dd), p_logrank = km$p)
  host_panels[[v]] <- km$panel
}
host_res <- bind_rows(host_res)
write.csv(host_res, "3-Output/04_host_variables_logrank.csv", row.names = FALSE)
cat("\nHost variables, log-rank PFS:\n"); print(host_res %>% mutate(p_logrank = signif(p_logrank, 3)), row.names = FALSE)
ggsave("2-Figure/FigS4_PFS_by_host_variables.svg",
       ggarrange(plotlist = host_panels, nrow = 5, ncol = 5), width = 18, height = 22)
cat("\nDone 04\n")
