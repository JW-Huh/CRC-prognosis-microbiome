################################################################################
##  Gut microbiota and CRC prognosis  (Huh et al., Microbiome 2022;10:203)
##  00_setup.R  -  shared environment for scripts 01-06
################################################################################
##  - Sourced at the top of every analysis script; not meant to be run alone
##  - Working directory = repository root (the folder holding 1-Input/, 4-Script/)
##  - Cohort : 333 primary CRC patients (SNUH), stool within 2 weeks before surgery
##  - Data   : 16S rRNA V3-V4 (MiSeq 2x300), QIIME2 + DADA2, SILVA 138 (99 % OTU
##             reference trained on V3-V4); see README for the upstream pipeline
##  - Outcome: PFS (progression or death), OS; follow-up ~ 2 years (PFS) / 3.5 y (OS)
################################################################################

## 0. Working directory and output folders -------------------------------------
##  - every path in this repository is relative to the repository root
if (!all(dir.exists(c("1-Input", "4-Script"))))
  stop("Set the working directory to the repository root ",
       "(the folder that contains 1-Input/ and 4-Script/).")
invisible(lapply(c("2-Figure", "3-Output"), dir.create, showWarnings = FALSE))

## 1. Packages ------------------------------------------------------------------
##  - tidyverse / magrittr   : data handling, %>%, set_colnames()
##  - vegan                  : rarefaction, alpha diversity, Bray-Curtis, PERMANOVA
##  - survival / survminer   : log-rank test, Cox PH, Kaplan-Meier figures
##  - ggpubr / ggrepel / ggforce : panel assembly, text labels, sina plots
suppressPackageStartupMessages({
  library(tidyverse); library(magrittr); library(vegan)
  library(survival);  library(survminer)
  library(ggpubr);    library(ggrepel);  library(ggforce)
})
options(stringsAsFactors = FALSE, dplyr.summarise.inform = FALSE, warn = 1)

## 2. Constants used across scripts --------------------------------------------
RAREFY_DEPTH <- 64604   # minimum library size among 333 samples = rarefaction depth
N_SAMPLE     <- 333
PREV_MIN     <- 10      # % prevalence : taxon present in >= 10 % of samples
ABUND_MIN    <- 0.01    # % abundance  : taxon reaches >= 0.01 % of all rarefied reads
SEED_RAREFY  <- 123     # vegan::rrarefy()
SEED_KMEANS  <- 12345   # k-means enterotype
SEED_CV      <- 123     # Monte-Carlo cross-validation and bootstrap of hazard scores
SEED_PERM    <- 1       # PERMANOVA (adonis2)

## 3. Colour codes ---------------------------------------------------------------
##  - enterotype 1 = Bacteroides-type (red), enterotype 2 = Prevotella-type (green)
col_enterotype <- c("1" = "red3", "2" = "green4")
col_stage      <- c("0" = "cornsilk1", "1" = "cornsilk2", "2" = "cornsilk3",
                    "3" = "coral1",    "4" = "coral4")
col_stage_cat  <- c("Developing" = "cornsilk2", "Advanced" = "coral4")
col_YN         <- c("N" = "cornsilk2", "Y" = "coral4")
col_01         <- c("0" = "cornsilk2", "1" = "coral4")
col_recur      <- c("yet" = "cornsilk2", ">1y" = "coral2", "<1y" = "coral4")
col_family     <- c("Prevotellaceae" = "green4",   "Bacteroidaceae"  = "red3",
                    "Ruminococcaceae" = "#D39200", "Lachnospiraceae" = "bisque",
                    "Oscillospiraceae" = "#93AA00", "Enterobacteriaceae" = "#F8766D",
                    "Akkermansiaceae" = "coral4",  "Fusobacteriaceae" = "#DB72FB",
                    "Rikenellaceae" = "steelblue", "Other" = "gray90")

## 4. Helper functions -----------------------------------------------------------

## 4-1. Relative abundance (%) at the rarefied depth
rel_pct <- function(x) x / RAREFY_DEPTH * 100

## 4-2. Dichotomise at the median
##  - "High" if strictly greater than the median, otherwise "Low"
##  - for sparse species (median = 0) this is simply detected vs. not detected
hi_lo <- function(x) ifelse(x > median(x, na.rm = TRUE), "High", "Low")

## 4-3. Clean taxon names so that they can be used as column names / formula terms
##  - "[Eubacterium]_eligens_group" -> "Eubacterium_eligens_group"
##  - "Escherichia-Shigella"        -> "Escherichia_Shigella"
clean_taxon <- function(x) x %>% gsub("\\[|\\]", "", .) %>% gsub("-", "_", .)

## 4-4. Log-rank test p-value for a survival formula
logrank_p <- function(formula, data) broom::glance(survdiff(formula, data = data))$p.value

## 4-5. One row of a Cox model (term = row index in the coefficient table)
cox_row <- function(formula, data, term = 1) {
  s <- summary(coxph(formula, data = data))
  data.frame(HR      = s$conf.int[term, 1],
             lower95 = s$conf.int[term, 3],
             upper95 = s$conf.int[term, 4],
             p       = s$coefficients[term, 5],
             C       = unname(s$concordance[1]),
             C_se    = unname(s$concordance[2]))
}

## 4-6. Kaplan-Meier panel in the style used throughout the paper
##  - curve (top) + number at risk (bottom), 5:1 height ratio
##  - axis text is removed on purpose; axis labels were typeset in Illustrator
##  - returns the survfit, the assembled panel and the log-rank p
km_plot <- function(formula, data, legend_labs, palette, file = NULL,
                    width = 5, height = 5.5, heights = c(5, 1),
                    legend_pos = c(0.88, 0.14)) {
  fit <- surv_fit(formula, data = data)
  g <- ggsurvplot(fit, data = data,
                  pval = TRUE, pval.coord = c(0, 0.05), size = 0.5,
                  censor.shape = "|", censor.size = 1,
                  legend.title = "", legend.labs = legend_labs, palette = palette,
                  risk.table = TRUE, risk.table.col = "strata", risk.table.y.text = FALSE)
  g$plot <- g$plot + theme_classic() +
    theme(axis.text = element_blank(), axis.title = element_blank(),
          axis.line = element_line(linewidth = 0.25),
          legend.text = element_text(size = rel(1.6)), legend.title = element_blank(),
          legend.background = element_rect(fill = NA),
          legend.position = "inside", legend.position.inside = legend_pos)
  g$table <- g$table + theme_classic() +
    theme(title = element_blank(), legend.position = "none",
          axis.line = element_line(linewidth = 0.25), axis.text.y = element_blank(),
          axis.text.x = element_text(size = rel(1.2), colour = "black"))
  panel <- ggarrange(g$plot, g$table, nrow = 2, ncol = 1, heights = heights)
  if (!is.null(file)) ggsave(file, panel, width = width, height = height)
  invisible(list(fit = fit, panel = panel, p = surv_pvalue(fit, data = data)$pval))
}

## 4-7. PCoA scatter (PCo1 vs PCo2) coloured by a variable, PERMANOVA p annotated
##  - my_var is passed unquoted (tidy evaluation); my_p is a number or NA
draw_pcoa <- function(my_df, my_var, my_p = NA, p_xy = c(-0.33, 0.24)) {
  p <- my_df %>%
    filter(!is.na({{ my_var }})) %>%
    ggplot(aes(PCo1, PCo2, fill = {{ my_var }})) +
    geom_point(shape = 21, size = 3, alpha = 0.9, stroke = 0.5) +
    theme_classic() +
    theme(legend.position = "top", axis.title = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(), aspect.ratio = 1)
  if (!is.na(my_p))
    p <- p +
      annotate("text", p_xy[1], p_xy[2], size = 6, label = "paste(italic(p), \"=\")", parse = TRUE) +
      annotate("text", p_xy[1] + 0.085, p_xy[2] + 0.006, size = 6, label = signif(my_p, 2))
  p
}

## 4-8. PERMANOVA p-value (Bray-Curtis, 999 permutations) for one grouping variable
##  - samples with NA in the variable are dropped before the test
permanova_p <- function(dist_obj, data, var, seed = SEED_PERM) {
  keep <- !is.na(data[[var]])
  ids  <- as.character(data$NO.[keep])
  dd   <- as.dist(as.matrix(dist_obj)[ids, ids])
  f    <- as.formula(paste0("dd ~ `", var, "`"))
  set.seed(seed)
  adonis2(f, data = data[keep, ])$`Pr(>F)`[1]
}

## 4-9. Waffle plot (one tile = one patient); replaces the 'waffle' package
##  - counts : named integer vector, category -> n; tiles fill column-wise
waffle_plot <- function(counts, colours, rows = 5) {
  cat_vec <- rep(names(counts), counts)
  df <- data.frame(cat = factor(cat_vec, levels = names(counts)),
                   idx = seq_along(cat_vec) - 1) %>%
    mutate(x = idx %/% rows, y = idx %% rows)
  ggplot(df, aes(x, y, fill = cat)) +
    geom_tile(colour = "white", linewidth = 0.8) +
    scale_fill_manual(values = colours) +
    coord_equal() + theme_void() + theme(legend.position = "none")
}

## 4-10. Upper triangle of a square matrix in long form (Source, Target, Value)
flatten_upper <- function(mat) {
  ut <- upper.tri(mat)
  data.frame(Source = rownames(mat)[row(mat)[ut]],
             Target = rownames(mat)[col(mat)[ut]],
             Value  = mat[ut])
}

## 5. Minimal ggplot themes -----------------------------------------------------
t_none  <- theme(legend.position = "none")
t_clean <- theme_classic() + theme(legend.position = "none")
t_nox   <- theme_classic() + theme(legend.position = "none", axis.title.x = element_blank())
