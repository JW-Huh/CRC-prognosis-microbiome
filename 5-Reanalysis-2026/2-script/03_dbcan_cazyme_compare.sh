#!/usr/bin/env bash
################################################################################
##  03_dbcan_cazyme_compare.sh   (TEMPLATE - drafted Oct 2026, not yet run)
##  CAZyme repertoire of Segatella copri vs Leyella stercorea (dbCAN3)
################################################################################
##  Question: does Leyella carry pectin / mannan loci that Segatella lacks? That
##  would give a genomic basis for the Leyella - boiled stem vegetable (Food46) link.
##  Environment: Linux / WSL with conda (needs internet for genomes + dbCAN DB)
##  Time: DB build ~20 min, annotation ~5 min per genome
################################################################################
set -euo pipefail
WD="$HOME/dbcan_prevotellaceae"; mkdir -p "$WD"/{genomes,db,out}; cd "$WD"

## 0. one-time install
# conda create -n dbcan -c bioconda -c conda-forge dbcan prodigal ncbi-datasets-cli -y
conda activate dbcan

## 1. reference genomes (RefSeq representative assemblies as of Oct 2026)
##    Segatella copri DSM 18205 (type), Leyella stercorea DSM 18206 (type, ex P. stercorea CB35)
##    Bacteroides thetaiotaomicron VPI-5482 as a broad-repertoire control
##    If an accession changed:  datasets summary genome taxon "Segatella copri" --reference
declare -A ACC=(
  [Segatella_copri]=GCF_020735445.1
  [Leyella_stercorea]=GCF_000235885.1
  [Bacteroides_theta]=GCF_014131755.1
)
for k in "${!ACC[@]}"; do
  datasets download genome accession "${ACC[$k]}" --include protein,genome --filename "genomes/$k.zip"
  unzip -o "genomes/$k.zip" -d "genomes/$k"
done

## 2. dbCAN database (one-time)
dbcan_build --cpus 8 --db-dir db --clean

## 3. annotation: protein FASTA if available, otherwise predict genes with prodigal
for k in "${!ACC[@]}"; do
  FAA=$(find "genomes/$k" -name "*.faa" | head -1)
  if [ -z "$FAA" ]; then
    FNA=$(find "genomes/$k" -name "*.fna" | head -1)
    prodigal -i "$FNA" -a "out/$k.faa" -p single -q; FAA="out/$k.faa"
  fi
  run_dbcan "$FAA" protein --out_dir "out/$k" --db_dir db --tools hmmer diamond dbcansub --hmm_cpu 8
done

## 4. GH / PL / CE / CBM family counts (families called by >= 2 tools)
python3 - <<'PY'
import pandas as pd, glob
rows = {}
for f in glob.glob("out/*/overview.txt"):
    name = f.split("/")[1]
    df = pd.read_csv(f, sep="\t"); df = df[df["#ofTools"] >= 2]
    rows[name] = df["HMMER"].str.extractall(r"(GH\d+|PL\d+|CE\d+|CBM\d+)")[0].value_counts()
tab = pd.DataFrame(rows).fillna(0).astype(int)
tab["Seg_minus_Ley"] = tab.get("Segatella_copri", 0) - tab.get("Leyella_stercorea", 0)
tab.sort_values("Seg_minus_Ley").to_csv("out/CAZyme_family_comparison.csv")
focus = ["GH10", "GH43", "GH51", "GH67", "GH115",   # xylan / arabinoxylan (expected Segatella-rich)
         "GH26", "GH113", "GH5",                     # mannan
         "GH28", "PL1", "PL9", "GH105",              # pectin (stem vegetables)
         "GH13", "GH77"]                             # starch
print(tab.reindex([f for f in focus if f in tab.index]).to_string())
PY
echo "done: $WD/out/CAZyme_family_comparison.csv"
## Reading the result: a large GH10/43/51 excess confirms Segatella as the xylan
## specialist; a relative Leyella excess in GH28 / PL (pectin) or GH26 (mannan)
## supports the stem-vegetable hypothesis.
