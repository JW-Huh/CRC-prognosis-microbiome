# 5-Reanalysis-2026 — Diet, gut microbiota and long-term outcome after CRC surgery

Follow-up to Huh et al., *Microbiome* 2022;10:203, handed over from Ji-Won Huh to
Da-Young Kang on 15 Oct 2026 (slides: `261015_handover_Kang.pptx`; status deck:
`Diet - CRC prognosis.pptx`). Everything in this folder builds on the main pipeline
in `4-Script/` (run `4-Script/run_all.R` first; it writes `3-Output/01_processed.rds`).

## Question
Two papers from this cohort cover diet → microbiota (Hoang et al., BMC Cancer 2022;
115 patients, 106-item FFQ, cross-sectional) and microbiota → outcome (Huh et al. 2022).
Neither joins the two. Does habitual diet explain the survival advantage of the
Prevotella-rich enterotype, or do diet and microbiota act independently?

## START HERE
1. Run `2-script/02_retaxonomy_prevotella.R` (needs `dada2` and the SILVA 138.2 /
   GTDB r220 reference files in `1-Input/db/`). Every Segatella / Leyella number
   below rests on a *provisional* assignment and must be confirmed first.
2. Then `2-script/01_diet_prognosis_reanalysis.R` (needs the 540-patient clinical +
   FFQ spreadsheet, which is **not** in the repository — see "Data not in git").

## How the provisional taxonomy was produced, and its limits
- The 1,004 ASVs labelled *Prevotella* by SILVA 138 were aligned (Needleman–Wunsch,
  `2-script/reassign.py`) to the 56 Prevotellaceae type-strain V3–V4 sequences in
  `1-Input/tree/Prevotella_16S_V3V4.fas`; best hit ≥ 97 % identity = species call,
  mapped to the genera of Hitch et al. 2022 (Syst Appl Microbiol 45:126354).
- No external database, no bootstrap confidence; 404 ASVs < 94 % identity were left
  as "Prevotellaceae_unresolved"; an ASV whose true species is absent from the 56
  references is forced onto its nearest neighbour.
- Therefore: direction is probably right, the values are provisional.

## Two data sets, mostly different patients — never merge on NO
- 540-patient clinical + FFQ spreadsheet: surgeries 2014–2021, FFQ n = 525.
- 333-patient microbiome cohort: stool collected Oct 2017 – Aug 2019.
- `NO` is not a shared key (merging on it gives 2 % age agreement).
- `1-input/id_link_verified.csv`: 31 pairs matched on sex, age, haemoglobin, WBC,
  BMI, stage and events; held-out CEA agrees in 97 %, platelets 94 %, albumin 94 %.
  Survival times differ by a constant ~78 days (later data freeze), ~29/31 genuine.
- Diet × microbiome analyses must use `1-Input/diet/Dietary_data_2022_115_residual_from_Image06.csv`
  (keyed to microbiome NO., already energy-adjusted by the residual method — do not
  re-adjust). The spreadsheet is for diet → outcome only.

## Folders
| Folder | Contents |
|---|---|
| `1-input` | provisional crosswalk (`prevotella_reassignment.csv`), genus tables before/after the split, `merged_333_newtax.csv` (metadata + CLR + re-derived enterotype), verified ID links |
| `2-script` | `01` diet → prognosis (R, template), `02` re-taxonomy (R, template), `03` dbCAN3 CAZyme comparison (bash, template); helpers `cox.py` (Cox + log-rank, numpy/scipy only), `rdata.py` (reads R `save.image()` files from Python), `reassign.py`, `make_panels2.py` (rebuilds P0–P7) |
| `3-figure` | P0–P7 as SVG + PNG, one panel per file |
| `4-output` | Cox tables, Food46 × 70-genus screen, PICRUSt2 contributions by genus |

## Headline numbers as they stand (provisional)
- Reads once called *Prevotella*: 81.7 % Segatella, 16.3 % Leyella, 1.3 % oral *Prevotella*.
- Prognosis (n = 333, stage-adjusted): Segatella OS HR 0.52; Leyella OS HR 0.60
  (per SD, p = 0.04); oral *Prevotella* points the other way. Re-derived enterotype
  (PAM on Jensen–Shannon) OS HR 0.50 (0.25–0.98).
- Food46 (bracken / sweet-potato stalk / taro stem) vs PFS in n = 114: log-rank
  p = 5.5e-4, HR 0.18 (0.06–0.52); ranked 1st of 209 dietary variables, FDR 0.12.
- Food46 correlates with Leyella (rho 0.29, FDR 0.16, 1st of 70 genera), not Segatella.
- PICRUSt2: thiamine biosynthesis 13.3 % Segatella vs 2.0 % Leyella; histidine
  degradation ~0 % Segatella vs 4.5 % Leyella — the two 2022 mechanisms split.
- Food46 keeps its effect after adjusting for either genus → not mediation (as tested).

## Data not in git
- `CRC data (540 patients)_clinical_diet_221220.xlsx` — contains patient names,
  hospital registration numbers and surgery/death dates. Keep it on the lab server;
  set the environment variable `CRC_540_XLSX` to its path before running script 01.
- Reference databases for script 02 (`1-Input/db/`, ~1 GB) — download from Zenodo.
