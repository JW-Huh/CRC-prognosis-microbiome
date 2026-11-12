################################################################################
##  03_microbial_hazard_score.R
##  31 candidate hazard scores -> bootstrap -> M5 -> comparison with clinical markers
################################################################################
##  Input  : 3-Output/01_processed.rds, 3-Output/02_screening.rds
##  Output : 3-Output/03_bootstrap_31_scores.csv, 03_cindex_indicators.csv,
##           03_cindex_combined_models.csv, 03_Table2_by_M5.csv, 03_risk.rds
##           2-Figure/Fig3*, FigS2d-e
##  Paper  : Fig. 3 A-G, Table 2, Figure S2 D-E
##  Checks : M5 Low/Moderate/High = 74/206/53 patients, progression 2/36/26
##           bootstrap M5: median log-rank p 2.1e-5, C 0.68 +- 0.03, mean HR 2.04
##           Stage + M5 (scaled) aHR 2.07, p 7.8e-9, C 0.78 (baseline stage 0.69)
##           Stage + M5 category: aHR moderate 6.53, high 19.54
################################################################################

source("4-Script/00_setup.R")
P <- readRDS("3-Output/01_processed.rds")
S <- readRDS("3-Output/02_screening.rds")
d <- S$d

## 1. Score components ----------------------------------------------------------
##  - Pv  : 1 point when Prevotella is LOW (<= median, i.e. < 16.1 % rel. abundance)
##  - Bsp, Din, Fnu, Ppi : 1 point when the pathobiont is detected
##    (Bacteroides sp. [= Alistipes], D. invisus, F. nucleatum, P. piscolens)
risk <- d %>%
  transmute(NO., PFS, PFS_time, OS, OS_time, Stage_1, CEA_1, `lymphatic invasion`, g__Prevotella,
            Pv  = as.integer(!(g__Prevotella > median(g__Prevotella))),
            Bsp = as.integer(s__Bacteroides_sp. > 0),
            Din = as.integer(s__Dialister_invisus > 0),
            Fnu = as.integer(s__Fusobacterium_nucleatum > 0),
            Ppi = as.integer(s__Pyramidobacter_piscolens > 0))
components <- c("Pv", "Bsp", "Din", "Fnu", "Ppi")

## 2. The 31 candidate scores = every non-empty subset of the five components ------
##  - score = number of unfavourable components present (0 .. k)
##  - ordered as in Fig. 3B: singles, Pv + 1, pairs without Pv, Pv + 2, Pv + 3,
##    triples without Pv, four pathobionts, all five (= M5)
subsets <- unlist(lapply(1:5, function(k) combn(components, k, simplify = FALSE)), recursive = FALSE)
names(subsets) <- sapply(subsets, paste, collapse = "+")
has_pv <- sapply(subsets, function(s) "Pv" %in% s); k_sub <- sapply(subsets, length)
score_order <- c(names(subsets)[k_sub == 1],
                 names(subsets)[k_sub == 2 &  has_pv], names(subsets)[k_sub == 2 & !has_pv],
                 names(subsets)[k_sub == 3 &  has_pv],
                 names(subsets)[k_sub == 4 &  has_pv], names(subsets)[k_sub == 3 & !has_pv],
                 names(subsets)[k_sub == 4 & !has_pv], names(subsets)[k_sub == 5])
stopifnot(length(score_order) == 31)
for (nm in names(subsets)) risk[[nm]] <- rowSums(risk[, subsets[[nm]], drop = FALSE])
risk <- risk %>% mutate(M5 = `Pv+Bsp+Din+Fnu+Ppi`,
                        M5_cat = factor(case_when(M5 >= 3 ~ "High", M5 >= 1 ~ "Moderate", TRUE ~ "Low"),
                                        levels = c("Low", "Moderate", "High")))
print(with(risk, table(PFS, M5_cat)))               # 2/74, 36/206, 26/53 progressed

## 3. Bootstrap of the 31 scores (Fig. 3B) -----------------------------------------
##  - 1,000 random half samples (166 / 333) without replacement
##  - per score: univariate Cox -> likelihood-ratio p, HR per point, C-index
##  - single Prevotella is re-dichotomised within each half sample (High vs Low,
##    so its HR is < 1); every other score uses the whole-cohort cut-off
##  - stability of a score = variance of log(HR) across the 1,000 subsamples
boot <- array(NA_real_, dim = c(1000, 31, 3), dimnames = list(NULL, score_order, c("p", "HR", "C")))
set.seed(SEED_CV)
for (i in 1:1000) {
  sub <- sample_n(risk, floor(nrow(risk) / 2), replace = FALSE) %>%
    mutate(Pv = g__Prevotella > median(g__Prevotella))   # within-sample High/Low
  for (sc in score_order) {
    cx <- summary(coxph(as.formula(paste0("Surv(PFS_time, PFS) ~ `", sc, "`")), data = sub))
    boot[i, sc, ] <- c(cx$logtest["pvalue"], cx$conf.int[1, 1], cx$concordance["C"])
  }
}
boot_summary <- data.frame(score = score_order,
                           n_components = k_sub[score_order],
                           median_p  = apply(boot[, , "p"], 2, median),
                           mean_HR   = apply(boot[, , "HR"], 2, mean),
                           min_HR    = apply(boot[, , "HR"], 2, min),
                           var_logHR = apply(log(boot[, , "HR"]), 2, var),
                           mean_C    = apply(boot[, , "C"], 2, mean),
                           sd_C      = apply(boot[, , "C"], 2, sd))
write.csv(boot_summary, "3-Output/03_bootstrap_31_scores.csv", row.names = FALSE)
cat("\nM5 in bootstrap:\n"); print(boot_summary[boot_summary$score == "Pv+Bsp+Din+Fnu+Ppi", ])

## Fig. 3B: -log10 p (top), C-index (middle), var(log HR) (bottom) + composition tile
long <- function(k) as.data.frame(boot[, , k]) %>% pivot_longer(everything(), names_to = "score") %>%
  mutate(score = factor(score, levels = score_order))
p_p <- long("p") %>% ggplot(aes(score, -log10(value))) +
  geom_hline(yintercept = -log10(0.05), linewidth = 1, linetype = "dotted", colour = "gray70") +
  geom_violin(width = 1, colour = "transparent", fill = "gray85", trim = TRUE) +
  geom_boxplot(aes(fill = score), width = 0.4, alpha = 0.5, outlier.shape = NA) +
  theme_classic() + theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
p_c <- long("C") %>% ggplot(aes(score, value)) +
  geom_hline(yintercept = 0.5, linewidth = 1, linetype = "dotted", colour = "gray") +
  geom_violin(width = 1, colour = "transparent", fill = "gray85", trim = TRUE) +
  geom_boxplot(aes(fill = score), width = 0.4, alpha = 0.5, outlier.shape = NA) +
  theme_classic() + theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
p_v <- boot_summary %>% mutate(score = factor(score, levels = score_order), is_M5 = score == "Pv+Bsp+Din+Fnu+Ppi") %>%
  ggplot(aes(score, var_logHR, group = 1)) + geom_line() + geom_point(aes(colour = is_M5)) +
  scale_colour_manual(values = c("gray50", "red3")) +
  theme_classic() + theme(legend.position = "none", axis.text = element_blank(), axis.title = element_blank())
tile <- expand_grid(score = score_order, component = components) %>%
  mutate(present = mapply(function(s, c) c %in% subsets[[s]], score, component),
         score = factor(score, levels = score_order), component = factor(component, levels = rev(components)))
p_t <- ggplot(tile, aes(score, component, fill = present)) + geom_tile(colour = "white") +
  scale_fill_manual(values = c("white", "orange")) + theme_classic() +
  theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
ggsave("2-Figure/Fig3b_bootstrap_31_scores.svg",
       ggarrange(p_p, p_c, p_v, p_t, ncol = 1, heights = c(1.4, 1.4, 0.8, 0.7)), width = 14, height = 7)

## 4. M5 survival curves (Fig. 3C-D; six-level score in Figure S2 D-E) -----------
pal6 <- c("green4", "yellow3", "yellow4", "red2", "red3", "red4")
risk$M5_f <- factor(risk$M5)
km_plot(Surv(PFS_time, PFS) ~ M5_f, risk, levels(risk$M5_f), pal6, "2-Figure/FigS2d_PFS_M5_six_levels.svg",
        heights = c(5, 1.5), height = 6.5, legend_pos = c(0.88, 0.18))
km_plot(Surv(OS_time, OS) ~ M5_f, risk, levels(risk$M5_f), pal6, "2-Figure/FigS2e_OS_M5_six_levels.svg",
        heights = c(5, 1.5), height = 6.5, legend_pos = c(0.88, 0.18))
km_plot(Surv(PFS_time, PFS) ~ M5_cat, risk, c("Low", "Moderate", "High"), c("green4", "yellow3", "red2"),
        "2-Figure/Fig3c_PFS_M5_category.svg", heights = c(5, 1.5), height = 5.8)
km_plot(Surv(OS_time, OS) ~ M5_cat, risk, c("Low", "Moderate", "High"), c("green4", "yellow3", "red2"),
        "2-Figure/Fig3d_OS_M5_category.svg", heights = c(5, 1.5), height = 5.8)

## 5. M5 versus established prognostic indicators (Fig. 3E-F) ---------------------
##  - baseline model: clinicopathological stage (0-2 vs 3/4)
##  - add one indicator at a time; report its adjusted HR and the model C-index
##  - continuous indicators are z-scaled (HR per SD); categorical ones use the
##    published cut-offs (High vs Low); enterotype = type 2 over type 1
risk2 <- d %>%
  transmute(PFS, PFS_time, Stage_1,
            Age, BMI, NLR, fibrinogen, CEA, g__Prevotella, M5 = risk$M5,
            Age_1, BMI_1, Sex, NLR_1, fibrinogen_1, CEA_1, `lymphatic invasion` = as.character(`lymphatic invasion`),
            Cluster_2, g__Prevotella_1, M5_cat = risk$M5_cat) %>%
  mutate(across(c(Age, BMI, NLR, fibrinogen, CEA, g__Prevotella, M5), ~ as.numeric(scale(.x))),
         across(c(Age_1, BMI_1, NLR_1, fibrinogen_1, CEA_1, g__Prevotella_1), ~ factor(.x, levels = c("Low", "High"))))
cx_base <- cox_row(Surv(PFS_time, PFS) ~ Stage_1, risk2)
cat("\nBaseline (stage) C-index:", round(cx_base$C, 3), "\n")

cx_num <- c("Age", "BMI", "NLR", "fibrinogen", "CEA", "g__Prevotella", "M5")
cx_cat <- c("Age_1", "BMI_1", "Sex", "NLR_1", "fibrinogen_1", "CEA_1", "lymphatic invasion",
            "Cluster_2", "g__Prevotella_1", "M5_cat")
label_of <- c(Age = "Age", BMI = "BMI", NLR = "NLR", fibrinogen = "Fibrinogen", CEA = "CEA",
              g__Prevotella = "Prevotella", M5 = "M5", Age_1 = "Age", BMI_1 = "BMI", Sex = "Sex",
              NLR_1 = "NLR", fibrinogen_1 = "Fibrinogen", CEA_1 = "CEA", `lymphatic invasion` = "Invasion",
              Cluster_2 = "Enterotype", g__Prevotella_1 = "Prevotella", M5_cat = "M5")
cindex <- lapply(c(cx_num, cx_cat), function(v) {
  f <- as.formula(paste0("Surv(PFS_time, PFS) ~ Stage_1 + `", v, "`"))
  nt <- if (v == "M5_cat") 2:3 else 2                 # M5 category has two contrasts
  lapply(nt, function(t) cox_row(f, risk2, term = t) %>%
           mutate(Variable = v, ID = label_of[v], Type = ifelse(v %in% cx_num, "Continuous", "Categorized"),
                  Level = if (v == "M5_cat") c("Moderate", "High")[t - 1] else "High vs Low")) %>% bind_rows()
}) %>% bind_rows() %>% select(Variable, ID, Type, Level, everything())
write.csv(cindex, "3-Output/03_cindex_indicators.csv", row.names = FALSE)
cat("\nStage + indicator:\n"); print(cindex %>% select(ID, Type, Level, HR, p, C) %>% mutate(across(c(HR, C), ~ round(.x, 3))))

cx_ID <- cindex %>% filter(Type == "Categorized") %>% arrange(-C) %>% pull(ID) %>% unique()
p_c <- cindex %>% mutate(ID = factor(ID, levels = rev(cx_ID))) %>%
  ggplot(aes(ID, C, colour = Type, fill = Type, group = Type)) +
  geom_hline(yintercept = cx_base$C, linewidth = 1, colour = "gray60", linetype = "dotted") +
  geom_errorbar(aes(ymin = C - C_se, ymax = C + C_se), position = "dodge", width = 0.5) +
  geom_point(aes(shape = Type), size = 4, colour = "black", position = position_dodge2(preserve = "single", width = 0.75)) +
  scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = c("gray40", "gray90")) +
  scale_colour_manual(values = c("gray40", "gray90")) + theme_classic() + coord_flip() +
  theme(legend.position = "none", axis.title.x = element_text(size = rel(1.4)),
        axis.text.y = element_text(size = rel(2), colour = "black"), axis.text.x = element_text(size = rel(1.1))) +
  labs(x = "", y = "Concordance index")
ggsave("2-Figure/Fig3e_Cindex_indicators.svg", p_c, width = 6, height = 4.5)

cx_ID <- cindex %>% filter(Type == "Categorized") %>% arrange(-HR) %>% pull(ID) %>% unique()
p_h <- cindex %>% mutate(ID = factor(ID, levels = rev(cx_ID))) %>%
  ggplot(aes(ID, log2(HR), colour = Type, fill = Type, group = Type)) +
  geom_hline(yintercept = 0, linewidth = 1, colour = "gray60", linetype = "dotted") +
  geom_errorbar(aes(ymin = log2(lower95), ymax = log2(upper95)), position = "dodge", width = 0.5) +
  geom_point(aes(shape = Type), size = 4, colour = "black", position = position_dodge2(preserve = "single", width = 0.75)) +
  scale_shape_manual(values = c(21, 23)) + scale_fill_manual(values = c("gray40", "gray90")) +
  scale_colour_manual(values = c("gray40", "gray90")) + theme_classic() + coord_flip() +
  theme(legend.position = "none", axis.title.x = element_text(size = rel(1.4)),
        axis.text.y = element_text(size = rel(2), colour = "black"), axis.text.x = element_text(size = rel(1.1))) +
  labs(x = "", y = "Log2 HR")
ggsave("2-Figure/Fig3f_HR_indicators.svg", p_h, width = 6, height = 4.5)

## 6. Does M5 add to a multivariable clinical model? (Fig. 3G) --------------------
##  - 15 nested models built from stage, CEA, lymphatic invasion, age, M5
##  - delta C-index of adding M5 ~ adding CEA + invasion + age together (+0.088)
cx_var <- c("Stage", "CEA", "Lym. Inv.", "Age", "M5")
models <- list("Stage" = "Stage_1", "CEA" = "CEA_1", "Lym. Inv." = "`lymphatic invasion`", "Age" = "Age_1", "M5" = "M5_cat",
               "Stage + CEA" = "Stage_1 + CEA_1", "Stage + Lym. Inv." = "Stage_1 + `lymphatic invasion`",
               "Stage + Age" = "Stage_1 + Age_1", "Stage + M5" = "Stage_1 + M5_cat",
               "Stage + CEA + Lym. Inv." = "Stage_1 + CEA_1 + `lymphatic invasion`",
               "Stage + CEA + Age" = "Stage_1 + CEA_1 + Age_1", "Stage + CEA + M5" = "Stage_1 + CEA_1 + M5_cat",
               "Stage + CEA + Lym. Inv. + Age" = "Stage_1 + CEA_1 + `lymphatic invasion` + Age_1",
               "Stage + CEA + Lym. Inv. + M5" = "Stage_1 + CEA_1 + `lymphatic invasion` + M5_cat",
               "Stage + CEA + Lym. Inv. + Age + M5" = "Stage_1 + CEA_1 + `lymphatic invasion` + Age_1 + M5_cat")
cx_combo <- data.frame(model = names(models),
                       C = sapply(models, function(f) summary(coxph(as.formula(paste("Surv(PFS_time, PFS) ~", f)), data = risk2))$concordance[1]))
write.csv(cx_combo, "3-Output/03_cindex_combined_models.csv", row.names = FALSE)
cat("\nCombined models:\n"); print(cx_combo %>% mutate(C = round(C, 4)), row.names = FALSE)
cat("delta C  Stage+M5 vs Stage:", round(cx_combo$C[9] - cx_combo$C[1], 4),
    "; Stage+CEA+Inv+Age vs Stage:", round(cx_combo$C[13] - cx_combo$C[1], 4), "\n")

tile <- expand_grid(model = names(models), var = cx_var) %>%
  mutate(present = mapply(function(m, v) v %in% strsplit(m, " \\+ ")[[1]], model, var),
         model = factor(model, levels = rev(names(models))), var = factor(var, levels = cx_var))
p_t <- ggplot(tile, aes(var, model, fill = present)) + geom_tile(colour = "white") +
  scale_fill_manual(values = c("white", "orange")) + theme_classic() +
  theme(legend.position = "none", axis.text = element_blank(), axis.ticks = element_blank(), axis.title = element_blank())
p_d <- cx_combo %>% mutate(model = factor(model, levels = rev(names(models)))) %>%
  ggplot(aes(C, model)) + geom_point() + theme_classic() +
  theme(axis.title = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())
ggsave("2-Figure/Fig3g_Cindex_combined_models.svg", ggarrange(p_t, p_d, ncol = 2, widths = c(1, 1.5)), width = 6, height = 6)

## 7. Table 2: clinical characteristics by M5 category -----------------------------
##  - continuous: one-way ANOVA; categorical: Fisher's exact test, or chi-square
##    where Fisher is not feasible (ASA, tumour location, stage)
##  - drug columns are coded NA (not taking) / 1 (taking); MSI 1 = MSS, 2/3 = MSI-L/H
##  - Smoking / Alcohol: 0 = never, 1 = former ("was"), 2 = current (Table 1: 47 smokers, 18 former);
##    tested as three levels, as in the paper
##  - MSI (26 missing) and K-ras (35 missing): missing patients are left out of the test (n/N)
##  - NOTE Table 2 of the paper prints "30 (57.7 %)" CEA-high patients in the M5-high group;
##    the data give 22 of 52 (42.3 %) and the published P value (1.2e-4) corresponds to 22
t2 <- d %>% mutate(M5_cat = risk$M5_cat,
                   Smoking = factor(Smoking, levels = 0:2, labels = c("never", "former", "current")),
                   Drinking = factor(Alcohol, levels = 0:2, labels = c("never", "former", "current")),
                   CEA_high = CEA_1 == "High", Lym_inv = `lymphatic invasion` == 1,
                   MSI_pos = ifelse(is.na(MSI), NA, MSI %in% c(2, 3)), Kras_mut = ifelse(is.na(Kras), NA, Kras == 1),
                   Medication = preoperative_medication == "1",
                   anti_HTN = !is.na(`antihypertensive drug`), anti_DM = !is.na(`antidiabetic drug`),
                   NSAID = !is.na(`Nonsteroidal anti-inflammatory drug`),
                   Location = factor(Location, levels = 1:4, labels = c("Proximal", "Distal", "Rectum", "Multiple")))
summ_num <- function(v) t2 %>% group_by(M5_cat) %>%
  summarise(value = sprintf("%.1f +- %.1f (%s-%s)", mean(.data[[v]], na.rm = TRUE), sd(.data[[v]], na.rm = TRUE),
                            min(.data[[v]], na.rm = TRUE), max(.data[[v]], na.rm = TRUE))) %>%
  pivot_wider(names_from = M5_cat, values_from = value) %>%
  mutate(Variable = v, p = summary(aov(t2[[v]] ~ t2$M5_cat))[[1]][["Pr(>F)"]][1], test = "ANOVA")
summ_cat <- function(v, test = "fisher") {
  tab <- table(t2[[v]], t2$M5_cat)
  p <- if (test == "fisher") fisher.test(tab, workspace = 2e7)$p.value else chisq.test(tab)$p.value
  as.data.frame.matrix(tab) %>% rownames_to_column("Level") %>%
    mutate(across(c(Low, Moderate, High), ~ sprintf("%d (%.1f%%)", .x, .x / sum(.x) * 100))) %>%
    mutate(Variable = v, p = p, test = test)
}
table2 <- bind_rows(summ_num("Age"), summ_cat("Sex"), summ_num("BMI"),
                    summ_cat("Smoking"), summ_cat("Drinking"), summ_cat("ASA", "chisq"),
                    summ_cat("CEA_high"), summ_cat("Comorbidity"), summ_cat("Hypertension"), summ_cat("Diabetes"),
                    summ_cat("Pulmonary disease"), summ_cat("Medication"), summ_cat("anti_HTN"), summ_cat("anti_DM"),
                    summ_cat("NSAID"), summ_cat("Location", "chisq"), summ_cat("Stage", "chisq"),
                    summ_cat("Lym_inv"), summ_cat("MSI_pos"), summ_cat("Kras_mut")) %>%
  select(Variable, Level, Low, Moderate, High, p, test)
write.csv(table2, "3-Output/03_Table2_by_M5.csv", row.names = FALSE)
cat("\nTable 2 (by M5 category):\n"); print(table2 %>% mutate(p = signif(p, 3)), row.names = FALSE)

## 8. Save -------------------------------------------------------------------------
saveRDS(list(risk = risk, boot = boot, boot_summary = boot_summary, cindex = cindex, cx_combo = cx_combo),
        "3-Output/03_risk.rds")
cat("\nSaved 3-Output/03_risk.rds\n")
