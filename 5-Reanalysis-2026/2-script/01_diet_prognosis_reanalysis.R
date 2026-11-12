################################################################################
##  01_diet_prognosis_reanalysis.R   (TEMPLATE - drafted Oct 2026, not yet run)
##  Diet -> Prevotella (Segatella / Leyella) -> CRC prognosis
################################################################################
##  Written for the 2026 follow-up (handover to Da-Young Kang, 15 Oct 2026).
##  Run from the repository root after the main pipeline (4-Script/01 ... 06).
##
##  Inputs
##   - <local>/CRC data (540 patients)_clinical_diet_221220.xlsx
##       clinical + FFQ for 540 surgical patients. NOT in the repository (contains
##       names, hospital IDs and dates). Point XL_540 at your local copy.
##   - 1-Input/diet/Dietary_variable_codebook.csv       variable names (EN / KR / code)
##   - 1-Input/diet/Dietary_data_2022_115_residual_from_Image06.csv
##       the 115-patient FFQ used in 2022, already energy-adjusted (residual method)
##   - 5-Reanalysis-2026/1-input/genus_newtax_counts.csv   genus counts after the
##       provisional Prevotella split (Segatella / Leyella / ...), rows = NO.
##   - 5-Reanalysis-2026/1-input/genus_oldtax_counts.csv   genus counts, SILVA 138 (2022)
##   - 3-Output/01_processed.rds                        master table d (333 patients)
##
##  Flow
##   Step 0  compare the 2022 FFQ file with the 540-patient spreadsheet (the two
##           are NOT the same patients; merging on NO is wrong - see README)
##   Step 1  FFQ pre-processing: residual-method energy adjustment (log scale)
##   Step 2  dietary patterns (17 food groups -> PCA, varimax)
##   Step 3  Cox: Food46 model table + dietary patterns -> PFS / OS  (n ~ 519)
##   Step 4  CLR Segatella / Leyella vs diet                         (n ~ 114)
##   Step 5  Prevotella HR with and without diet adjustment, interaction
##   Step 6  figures
################################################################################

rm(list = ls())
suppressPackageStartupMessages({
  library(tidyverse); library(readxl); library(survival); library(survminer)
  library(psych); library(ggrepel); library(patchwork)
})
if (!dir.exists("4-Script")) stop("run from the repository root")
dir.create("5-Reanalysis-2026/2-figure", showWarnings = FALSE)
dir.create("5-Reanalysis-2026/4-output", showWarnings = FALSE)
XL_540 <- Sys.getenv("CRC_540_XLSX", "../../CRC data (540 patients)_clinical_diet_221220.xlsx")

## Step 0. 2022 FFQ file vs the 540-patient spreadsheet --------------------------
##  - the 540 spreadsheet (surgeries 2014-2021) and the 333 microbiome patients
##    (stool Oct 2017 - Aug 2019) overlap only partly; NO is not a shared key
##  - 1-input/id_link_verified.csv holds 31 pairs matched on sex / age / labs
xl <- read_xlsx(XL_540) %>%
  mutate(across(c(age, sex, BMI, AJCCstage, OS, OS_time, PFS, PFS_time, DFS, DFS_time,
                  starts_with("food"), starts_with("var"), starts_with("group"), wet),
                ~ suppressWarnings(as.numeric(.x))))
code <- read.csv("1-Input/diet/Dietary_variable_codebook.csv", encoding = "UTF-8") %>%
  set_names(c("en", "ko", "var")) %>% filter(var != "")

old <- read.csv("1-Input/diet/Dietary_data_2022_115_residual_from_Image06.csv", check.names = FALSE)
cat("2022 FFQ file:", nrow(old), "rows,", ncol(old), "cols\n")
cat("2022 Food46 (energy-adjusted):\n"); print(summary(old$Food46))
cat("Spreadsheet food46 (raw g/day):\n"); print(summary(xl$food46))
link <- read.csv("5-Reanalysis-2026/1-input/id_link_verified.csv")
cmp  <- old %>% select(NO_mb = `NO.`, Food46_2022 = Food46) %>%
  inner_join(link, by = "NO_mb") %>%
  inner_join(xl %>% select(NO_xl = NO, food46), by = "NO_xl")
cat("verified pairs with both FFQ values:", nrow(cmp), "\n")
if (nrow(cmp) > 5) print(cor(cmp$Food46_2022, cmp$food46, method = "spearman", use = "complete"))
##  - a weak correlation here is expected: the 2022 values are energy residuals,
##    the spreadsheet values are raw grams; a near-zero correlation means the
##    two files describe different patients and must not be pooled

## Step 1. FFQ pre-processing and energy adjustment (540-patient spreadsheet) ------
##  - analysis set: FFQ present, PFS known, stage known
##  - residual method on the log scale: log(x + 1) ~ log(energy); residuals are
##    centred at zero, so add back the mean of log(x + 1) for interpretable units
##  - implausible energy (< 500 or > 4,000 kcal/day) excluded
dat <- xl %>%
  filter(!is.na(food46), !is.na(PFS), !is.na(PFS_time), PFS_time > 0, !is.na(AJCCstage),
         var1 >= 500, var1 <= 4000) %>%
  mutate(stage34  = as.integer(AJCCstage >= 3),
         surg_yr  = lubridate::year(as.Date(.data[[grep("surg|수술", names(xl), value = TRUE)[1]]])),
         early_yr = as.integer(surg_yr <= 2016),
         log_E    = log(var1))
cat("analysis n =", nrow(dat), " PFS events =", sum(dat$PFS), " OS events =", sum(dat$OS, na.rm = TRUE), "\n")

resid_adj <- function(x, logE) { y <- log1p(x); residuals(lm(y ~ logE)) + mean(y) }
grp_cols <- paste0("group", 1:17)
nut_cols <- c("var2", "var3", "var4", "var5", "var6", "var7")   # plant/animal protein, plant/animal fat, carbohydrate, fibre
adj <- dat %>% transmute(NO, across(all_of(c(grp_cols, nut_cols, "food46")), ~ resid_adj(.x, log_E), .names = "{.col}_r"))
dat <- left_join(dat, adj, by = "NO")

## Step 2. Dietary patterns (17 food groups, energy-adjusted -> PCA + varimax) ---------
##  - name the patterns from the loadings (e.g. "plant / traditional" vs "Western")
X   <- dat %>% select(all_of(paste0(grp_cols, "_r"))) %>% scale()
pca <- principal(X, nfactors = 3, rotate = "varimax", scores = TRUE)
print(pca$loadings, cutoff = 0.3)
dat$PC1 <- pca$scores[, 1]; dat$PC2 <- pca$scores[, 2]; dat$PC3 <- pca$scores[, 3]
write.csv(unclass(pca$loadings), "5-Reanalysis-2026/4-output/diet_pattern_loadings.csv")

## Step 3. Cox models ---------------------------------------------------------------
tidy_cox <- function(f, data, term = 1) {
  s <- summary(coxph(f, data = data))
  tibble(HR = s$conf.int[term, 1], lo = s$conf.int[term, 3], hi = s$conf.int[term, 4], p = s$coefficients[term, 5])
}
hi_lo <- function(x) as.integer(x >= median(x, na.rm = TRUE))
dat <- dat %>% mutate(f46_hi = hi_lo(food46), f46r_hi = hi_lo(food46_r),
                      f46_tert = ntile(food46, 3), f46_log = log1p(food46))

## 3a. Food46 under increasingly adjusted models (the design-document table)
models <- list(
  "crude"                          = Surv(PFS_time, PFS) ~ f46_hi,
  "+stage (as in 2022)"            = Surv(PFS_time, PFS) ~ f46_hi + stage34,
  "+stage +logE"                   = Surv(PFS_time, PFS) ~ f46_hi + stage34 + log_E,
  "+stage +logE +age +sex"         = Surv(PFS_time, PFS) ~ f46_hi + stage34 + log_E + age + sex,
  "residual hi/lo +stage +age +sex" = Surv(PFS_time, PFS) ~ f46r_hi + stage34 + age + sex,
  "continuous log +all"            = Surv(PFS_time, PFS) ~ f46_log + stage34 + log_E + age + sex,
  "per tertile +all"               = Surv(PFS_time, PFS) ~ f46_tert + stage34 + log_E + age + sex,
  "+all +surgery year"             = Surv(PFS_time, PFS) ~ f46_hi + stage34 + log_E + age + sex + surg_yr)
tab_f46 <- bind_rows(
  imap_dfr(models, ~ tidy_cox(.x, dat) %>% mutate(model = .y, outcome = "PFS")),
  imap_dfr(models, ~ tidy_cox(update(.x, Surv(OS_time, OS) ~ .), dat) %>% mutate(model = .y, outcome = "OS")))
print(tab_f46, n = 20)
write.csv(tab_f46, "5-Reanalysis-2026/4-output/cox_food46_models.csv", row.names = FALSE)

## 3b. Dietary patterns -> PFS / OS (tertiles); sensitivity: surgery 2014-16, stage III only
pat_cox <- function(pc, data, outcome = "PFS") {
  data$tert <- ntile(data[[pc]], 3)
  f <- if (outcome == "PFS") Surv(PFS_time, PFS) ~ tert + stage34 + age + sex + surg_yr
       else                  Surv(OS_time,  OS)  ~ tert + stage34 + age + sex + surg_yr
  tidy_cox(f, data) %>% mutate(pattern = pc, outcome = outcome)
}
tab_pat <- bind_rows(
  map_dfr(c("PC1", "PC2", "PC3"), ~ bind_rows(pat_cox(.x, dat, "PFS"), pat_cox(.x, dat, "OS"))) %>% mutate(set = "all"),
  map_dfr(c("PC1", "PC2", "PC3"), ~ pat_cox(.x, filter(dat, early_yr == 1), "PFS")) %>% mutate(set = "surgery 2014-16"),
  map_dfr(c("PC1", "PC2", "PC3"), ~ pat_cox(.x, filter(dat, AJCCstage == 3), "PFS")) %>% mutate(set = "stage III"))
print(tab_pat); write.csv(tab_pat, "5-Reanalysis-2026/4-output/cox_diet_patterns.csv", row.names = FALSE)

## 3c. Energy-adjusted nutrients (tertiles) - supplementary
tab_nut <- map_dfr(paste0(nut_cols, "_r"), ~ pat_cox(.x, dat, "PFS")) %>%
  left_join(code %>% transmute(pattern = paste0(var, "_r"), en), by = "pattern")
print(tab_nut); write.csv(tab_nut, "5-Reanalysis-2026/4-output/cox_nutrients_resid.csv", row.names = FALSE)

## Step 4. Microbiome (re-assigned genera) vs diet, in the 2022 FFQ patients ----------
##  - this step uses the 115-patient 2022 file, because only those patients have
##    both FFQ and stool; the 540 spreadsheet is for diet -> outcome only
P   <- readRDS("3-Output/01_processed.rds"); d <- P$d
gen     <- read.csv("5-Reanalysis-2026/1-input/genus_newtax_counts.csv", row.names = 1, check.names = FALSE)
gen_old <- read.csv("5-Reanalysis-2026/1-input/genus_oldtax_counts.csv", row.names = 1, check.names = FALSE)
clr_tab <- function(m, pseudo = 0.5) { m <- as.matrix(m); m[m == 0] <- pseudo; t(apply(m, 1, function(r) log(r) - mean(log(r)))) }
clr_new <- clr_tab(gen); clr_old <- clr_tab(gen_old)
prev_cols <- grep("Prevotella|Segatella|Hoylesella|Leyella|Palleniella|Xylanibacter|Hallella", colnames(gen), value = TRUE)
print(prev_cols)

diet22 <- old %>% select(-NO_mb) %>% mutate(`NO.` = as.character(`NO.`))
mb <- as.data.frame(clr_new[, prev_cols, drop = FALSE]) %>% rownames_to_column("NO.") %>%
  mutate(Prevotella_old = clr_old[match(`NO.`, rownames(clr_old)), grep("Prevotella", colnames(clr_old))[1]]) %>%
  inner_join(diet22, by = "NO.") %>%
  inner_join(d %>% select(`NO.`, PFS, PFS_time, OS, OS_time, Stage_1, Age, Sex), by = "NO.") %>%
  mutate(stage34 = as.integer(Stage_1 == "Advanced"))
cat("diet x microbiome overlap n =", nrow(mb), "\n")
seg <- grep("Segatella", prev_cols, value = TRUE)[1]
cat("2022 Prevotella vs Segatella (Spearman):", cor(mb$Prevotella_old, mb[[seg]], method = "spearman"), "\n")

##  Food groups (already energy-adjusted in the 2022 file) and Food46 vs CLR taxa
diet_vars <- c(paste0("group", 1:17), "var7", "Food46")
assoc <- expand_grid(diet = diet_vars, taxon = c(prev_cols, "Prevotella_old")) %>%
  mutate(res = map2(diet, taxon, ~ broom::tidy(cor.test(mb[[.x]], mb[[.y]], method = "spearman", exact = FALSE)))) %>%
  unnest(res) %>% mutate(FDR = p.adjust(p.value, "fdr")) %>% arrange(p.value)
print(assoc, n = 30); write.csv(assoc, "5-Reanalysis-2026/4-output/diet_vs_prevotella_clr.csv", row.names = FALSE)

## Step 5. Segatella HR before / after diet adjustment, and interaction (exploratory) ---
##  - n ~ 114 with ~22 PFS events: treat as hypothesis-generating
mb$seg_hi <- hi_lo(mb[[seg]]); mb$f46_hi <- hi_lo(mb$Food46)
tab_med <- bind_rows(
  tidy_cox(Surv(PFS_time, PFS) ~ seg_hi + stage34, mb)                     %>% mutate(model = "Segatella + stage"),
  tidy_cox(Surv(PFS_time, PFS) ~ seg_hi + stage34 + f46_hi, mb)            %>% mutate(model = "+ Food46"),
  tidy_cox(Surv(PFS_time, PFS) ~ seg_hi + stage34 + var7, mb)              %>% mutate(model = "+ fibre"),
  tidy_cox(Surv(PFS_time, PFS) ~ seg_hi * f46_hi + stage34, mb, term = 4)  %>% mutate(model = "Segatella x Food46 (interaction)"))
print(tab_med); write.csv(tab_med, "5-Reanalysis-2026/4-output/cox_prevotella_diet_adjusted.csv", row.names = FALSE)

## Step 6. Figures ------------------------------------------------------------------
for (pc in c("PC1", "PC2", "PC3")) {
  dat$tert <- factor(ntile(dat[[pc]], 3), labels = c("T1", "T2", "T3"))
  p <- ggsurvplot(survfit(Surv(PFS_time, PFS) ~ tert, data = dat), data = dat, pval = TRUE, risk.table = TRUE,
                  palette = c("gray60", "steelblue", "coral3"), legend.title = pc)
  ggsave(sprintf("5-Reanalysis-2026/2-figure/Fig2_KM_%s.svg", pc),
         ggarrange(p$plot, p$table, nrow = 2, heights = c(5, 1.3)), width = 5, height = 6)
}
p3 <- mb %>% ggplot(aes(Food46, .data[[seg]])) + geom_point(alpha = .6) + geom_smooth(method = "lm", formula = y ~ x) +
  theme_classic() + labs(y = paste0("CLR ", seg), x = "Food46 (energy-adjusted)")
ggsave("5-Reanalysis-2026/2-figure/Fig3_Food46_vs_Segatella.svg", p3, width = 5, height = 4)
p4 <- tab_med %>% mutate(model = fct_rev(fct_inorder(model))) %>%
  ggplot(aes(HR, model)) + geom_vline(xintercept = 1, linetype = 2, colour = "gray50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = .2) + geom_point(size = 3) +
  scale_x_log10() + theme_classic() + labs(x = "HR for PFS (log scale)", y = NULL)
ggsave("5-Reanalysis-2026/2-figure/Fig4_forest_prevotella_adjusted.svg", p4, width = 6, height = 3)
cat("Done 01 (2026 reanalysis)\n")
