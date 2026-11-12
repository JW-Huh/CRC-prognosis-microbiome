################################################################################
##  06_diet_exploratory_2022.R
##  Food-frequency data (115 patients) vs prognosis and vs gut taxa - EXPLORATORY
################################################################################
##  Input  : 1-Input/diet/Dietary_data_2022_115_residual_from_Image06.csv
##           3-Output/01_processed.rds
##  Output : 3-Output/06_diet_vs_prognosis.csv, 06_Food46_vs_taxa_wilcoxon.csv,
##           06_Food46_species_network.csv ; 2-Figure/Diet-00*
##  Status : NOT part of the 2022 paper. Run in Dec 2022 after publication as a
##           first look at diet; it seeded the 2026 follow-up (5-Reanalysis-2026/)
##  Diet data provenance
##   - 106-item semi-quantitative FFQ collected by the SNUH / SNU Preventive Medicine
##     team (Hoang T, ..., Shin A. BMC Cancer 2022;22:656), 115 of the 333 patients
##   - 209 variables: 93 nutrients (var*), 17 food groups (group*), 99 food items (Food*)
##   - nutrients and food groups were already energy-adjusted by the residual method
##     before I received them (negative values occur in 157 / 209 columns);
##     var1 = total energy (kcal/day) is the raw intake  -> do NOT re-adjust
##   - keyed to the microbiome sample number NO. ; no patient identifiers
################################################################################

source("4-Script/00_setup.R")
P <- readRDS("3-Output/01_processed.rds"); d <- P$d; taxa_to_test <- P$taxa_to_test

## 1. Import and merge --------------------------------------------------------------
diet <- read.csv("1-Input/diet/Dietary_data_2022_115_residual_from_Image06.csv", check.names = FALSE)
stopifnot(all(diet$NO. == diet$NO_mb))
diet <- diet %>% select(-NO_mb) %>% mutate(NO. = as.character(NO.))
diet_vars <- setdiff(colnames(diet), "NO.")
cat("Diet table:", nrow(diet), "patients x", length(diet_vars), "variables\n")
cat("Columns with negative values:", sum(sapply(diet[diet_vars], function(x) any(x < 0, na.rm = TRUE))), "/", length(diet_vars), "\n")
cat("Energy (var1, kcal/day): "); print(summary(diet$var1))

dd <- merge(diet, d, by = "NO.")                     # 114 patients with FFQ + stool
cat("Patients with FFQ and microbiome:", nrow(dd), "; PFS events", sum(dd$PFS), "; deaths", sum(dd$OS), "\n")

## 2. Screen every dietary variable against PFS and OS --------------------------------
##  - intake dichotomised at the median; log-rank test + Cox adjusted for stage (0-2 vs 3/4)
##  - only 22 progression events: a hypothesis-generating screen, FDR is reported
diet_bin <- dd %>% mutate(across(all_of(diet_vars), ~ ifelse(.x >= median(.x, na.rm = TRUE), "1", "0")))
diet_res <- lapply(diet_vars, function(v) {
  tmp <- diet_bin %>% select(PFS, PFS_time, OS, OS_time, Stage_1, var = all_of(v))
  cbind(data.frame(var = v, PFS_p_logrank = logrank_p(Surv(PFS_time, PFS) ~ var, tmp)),
        cox_row(Surv(PFS_time, PFS) ~ var + Stage_1, tmp) %>% select(HR, lower95, upper95) %>% setNames(paste0("PFS_", c("HR", "lower95", "upper95"))),
        data.frame(OS_p_logrank = logrank_p(Surv(OS_time, OS) ~ var, tmp)),
        cox_row(Surv(OS_time, OS) ~ var + Stage_1, tmp) %>% select(HR, lower95, upper95) %>% setNames(paste0("OS_", c("HR", "lower95", "upper95"))))
}) %>% bind_rows() %>%
  mutate(PFS_FDR = p.adjust(PFS_p_logrank, "fdr"), OS_FDR = p.adjust(OS_p_logrank, "fdr")) %>%
  arrange(PFS_p_logrank)
write.csv(diet_res, "3-Output/06_diet_vs_prognosis.csv", row.names = FALSE)
cat("\nTop dietary variables for PFS:\n"); print(head(diet_res %>% select(var, PFS_p_logrank, PFS_FDR, PFS_HR, OS_p_logrank), 10), row.names = FALSE)

volcano <- function(hr, p, file, ylab) {
  g <- diet_res %>% ggplot(aes(log(.data[[hr]]), -log10(.data[[p]]))) +
    geom_hline(yintercept = 2, linetype = "dotted", colour = "gray50") + geom_vline(xintercept = 0, linetype = "dotted", colour = "gray50") +
    geom_point(shape = 21, size = 4, aes(fill = log2(.data[[hr]]))) + geom_text_repel(aes(label = var)) +
    annotate("text", x = -1.5, y = 2.12, label = "P = 0.01", colour = "gray50", size = 6) +
    scale_fill_gradient2(low = "steelblue", mid = "white", high = "coral3") + theme_classic() +
    theme(legend.position = "none", axis.text = element_text(size = rel(1.05), colour = "black"), axis.title = element_text(size = rel(1.4))) +
    labs(x = "Log (hazard ratio, stage-adjusted)", y = ylab)
  ggsave(file, g, width = 6, height = 6)
}
volcano("PFS_HR", "PFS_p_logrank", "2-Figure/Diet-001_logrank_PFS_volcano.svg", "-Log10 (log-rank P), PFS")
volcano("OS_HR",  "OS_p_logrank",  "2-Figure/Diet-002_logrank_OS_volcano.svg",  "-Log10 (log-rank P), OS")

## 3. Food46 = bracken / sweet-potato stalk / taro stem (boiled stem vegetables) --------
##  - ranked first for PFS (log-rank p = 5.5e-4, FDR 0.12); high intake -> fewer events
km_plot(Surv(PFS_time, PFS) ~ Food46, diet_bin, c("Low", "High"), c("black", "green3"),
        "2-Figure/Diet-003_Food46_PFS.svg", width = 4, height = 5.5)
print(with(diet_bin, table(Food46, PFS)))

## 4. Food46 (High / Low) vs. gut taxa ---------------------------------------------------
##  - Wilcoxon rank-sum per taxon (209 genera + species from script 01)
##  - keep taxa with median abundance >= 0.05 % in at least one group
f46 <- dd %>% mutate(Food46_bin = ifelse(Food46 >= median(Food46), "High", "Low"))
f46_res <- lapply(taxa_to_test, function(t) {
  x <- f46[[t]]; g <- f46$Food46_bin
  data.frame(Taxon = t, WilcoxP = wilcox.test(x ~ g)$p.value,
             Med_High = median(rel_pct(x[g == "High"])), Med_Low = median(rel_pct(x[g == "Low"])))
}) %>% bind_rows() %>% filter(Med_High >= 0.05 | Med_Low >= 0.05) %>%
  mutate(FDR = p.adjust(WilcoxP, "fdr"), Direction = ifelse(Med_High > Med_Low, "Increased", "Decreased")) %>%
  arrange(WilcoxP)
write.csv(f46_res, "3-Output/06_Food46_vs_taxa_wilcoxon.csv", row.names = FALSE)
cat("\nTaxa differing by Food46 intake (p < 0.05):\n"); print(f46_res %>% filter(WilcoxP < 0.05), row.names = FALSE)

p <- f46 %>% select(Food46_bin, g__Prevotella, g__Parabacteroides, g__Dialister) %>%
  pivot_longer(-Food46_bin, names_to = "Genus", values_to = "Count") %>% mutate(Genus = gsub("g__", "", Genus)) %>%
  ggplot(aes(Food46_bin, log10(rel_pct(Count) + 1))) +
  geom_jitter(aes(colour = Food46_bin), size = 2, width = 0.2) + geom_boxplot(aes(colour = Food46_bin), alpha = 0.7, outlier.alpha = 0) +
  facet_wrap(~ Genus, nrow = 1, scales = "free_y") + stat_compare_means(comparisons = list(c("High", "Low")), method = "wilcox") +
  scale_colour_manual(values = c(High = "green3", Low = "gray40")) + theme_classic() +
  theme(legend.position = "none", axis.text.x = element_text(size = rel(1.5), colour = "black"), strip.text = element_text(size = rel(1.1))) +
  labs(x = "Food46 intake", y = "log10 (relative abundance % + 1)")
ggsave("2-Figure/Diet-004_Food46_three_genera.svg", p, width = 6, height = 5)

## 5. Co-abundance network around Food46 (for Cytoscape) --------------------------------
##  - Spearman correlation among Food46 and all tested taxa
##  - species positively correlated with Food46 (FDR < 0.2) and their positive
##    partners (FDR < 0.2) form the exported edge list
cor_mat <- f46 %>% select(Food46, all_of(taxa_to_test)) %>% as.matrix() %>% Hmisc::rcorr(type = "spearman")
f46_cor <- cbind(flatten_upper(cor_mat$r), P = flatten_upper(cor_mat$P)$Value) %>%
  setNames(c("Source", "Target", "Spearman", "Pval")) %>% mutate(FDR = p.adjust(Pval, "fdr")) %>% arrange(FDR)
f46_bugs <- f46_cor %>% filter(Source == "Food46", grepl("s__", Target), Spearman > 0, FDR < 0.2) %>% pull(Target)
f46_2nd  <- union(f46_cor %>% filter(Source %in% c("Food46", f46_bugs), FDR < 0.2, grepl("s__", Target)) %>% pull(Target),
                  f46_cor %>% filter(Target %in% c("Food46", f46_bugs), FDR < 0.2, grepl("s__", Source)) %>% pull(Source))
f46_network <- f46_cor %>% filter(Source %in% f46_2nd, Target %in% f46_2nd, FDR < 0.2, Spearman > 0)
write.csv(f46_network, "3-Output/06_Food46_species_network.csv", row.names = FALSE, quote = FALSE)
cat("\nSpecies directly correlated with Food46 (FDR < 0.2):", f46_bugs, "\n")
cat("Network edges exported:", nrow(f46_network), "\n")
cat("\nDone 06\n")
