"""Provisional re-assignment of the 1,004 'Prevotella' ASVs (SILVA 138) to 56 Prevotellaceae
type-strain V3-V4 references, then to the genera of Hitch et al. 2022 (Segatella, Hoylesella,
Leyella, Xylanibacter, Hallella, Prevotella sensu stricto).

Usage (from the repository root):
    python 5-Reanalysis-2026/2-script/reassign.py 5-Reanalysis-2026/1-input/prevotella_reassignment.csv

- k-mer (8-mer Jaccard) prefilter -> 3 candidates per ASV -> global Needleman-Wunsch identity
- identity >= 0.97 species call, 0.94-0.97 '_lowconf', < 0.94 'Prevotellaceae_unresolved'
- NOT a substitute for DADA2 assignTaxonomy against SILVA 138.2 / GTDB (see README); ~9 min
"""
import pandas as pd, numpy as np, re, sys
U = sys.argv[2] if len(sys.argv) > 2 else "1-Input/"        # folder holding the 2022 inputs
OUT = sys.argv[1]

def read_fasta(p):
    d, k = {}, None
    for line in open(p):
        line = line.strip()
        if not line: continue
        if line.startswith(">"): k = line[1:]; d[k] = ""
        else: d[k] += line.upper().replace("-", "").replace("U", "T")
    return d

refs = read_fasta(U + "tree/Prevotella_16S_V3V4.fas")
asvseq = read_fasta(U + "asv_seq.fasta")
tax = pd.read_csv(U + "CRC_prognosis_Taxonomy.csv").set_index("Feature ID")
asv = pd.read_csv(U + "CRC_prognosis_ASV_table.csv", index_col=0)
genus = tax.Taxon.str.extract(r"g__([^;]+)")[0].str.strip()
prev_ids = genus[genus == "Prevotella"].index
reads = asv.sum(axis=1)
print("Prevotella ASVs:", len(prev_ids), "with seq:", sum(i in asvseq for i in prev_ids))

# ---------- Hitch et al. 2022 mapping (LPSN, 2026-10) ----------
NEW = {
 "Segatella":   "albensis baroniae bryantii buccae cerevisiae copri hominis maculosa oris oulorum paludivivens salivae".split(),
 "Hoylesella":  "buccalis enoeca loescheii marshii nanceiensis oralis pleuritidis saccharolytica shahii timonensis".split(),
 "Leyella":     "stercorea lascolaii".split(),
 "Xylanibacter":"brevis oryzae rara ruminicola".split(),
 "Hallella":    "bergensis colorans mizrahii multisaccharivorax seregens".split(),
}
sp2genus = {sp: g for g, l in NEW.items() for sp in l}
def ref_species(h):
    h = h.replace("Hallela ", "Prevotella_")
    return h.split("_")[1]
def ref_genus(h): return sp2genus.get(ref_species(h), "Prevotella")

# ---------- k-mer prefilter ----------
K = 8
def kmers(s): return set(s[i:i+K] for i in range(len(s)-K+1))
ref_names = list(refs); ref_k = [kmers(refs[r]) for r in ref_names]

# ---------- Needleman-Wunsch identity (vectorized per row) ----------
def nw_identity(a, b, match=1, mismatch=-1, gap=-2):
    n, m = len(a), len(b)
    A = np.frombuffer(a.encode(), dtype=np.uint8); B = np.frombuffer(b.encode(), dtype=np.uint8)
    S = np.where(A[:, None] == B[None, :], match, mismatch)
    H = np.zeros((n+1, m+1), dtype=np.int32); H[0, :] = gap*np.arange(m+1); H[:, 0] = gap*np.arange(n+1)
    P = np.zeros((n+1, m+1), dtype=np.int8)  # 0 diag,1 up,2 left
    for i in range(1, n+1):
        diag = H[i-1, :-1] + S[i-1]
        up = H[i-1, 1:] + gap
        row = np.maximum(diag, up); ptr = np.where(diag >= up, 0, 1)
        # left dependency: sequential within row
        left = H[i, 0]
        for j in range(1, m+1):
            l = left + gap
            if l > row[j-1]: row[j-1] = l; ptr[j-1] = 2
            left = row[j-1]
        H[i, 1:] = row; P[i, 1:] = ptr
    # traceback
    i, j, ident, aln = n, m, 0, 0
    while i > 0 or j > 0:
        if i > 0 and j > 0 and P[i, j] == 0:
            ident += (a[i-1] == b[j-1]); i -= 1; j -= 1
        elif i > 0 and (j == 0 or P[i, j] == 1): i -= 1
        else: j -= 1
        aln += 1
    return ident / aln

rows = []
for n, aid in enumerate(prev_ids):
    s = asvseq.get(aid)
    if s is None: continue
    ks = kmers(s)
    jac = np.array([len(ks & rk) / len(ks | rk) for rk in ref_k])
    top = np.argsort(-jac)[:3]
    best = max(((nw_identity(s, refs[ref_names[t]]), t) for t in top), key=lambda x: x[0])
    second = sorted(((nw_identity(s, refs[ref_names[t]]), t) for t in top), reverse=True)
    ident, t = best
    r = ref_names[t]
    rows.append(dict(ASV=aid, reads=int(reads.get(aid, 0)), best_ref=r, species=ref_species(r),
                     new_genus=ref_genus(r), identity=round(ident, 4),
                     second_ref=ref_names[second[1][1]] if len(second) > 1 else None,
                     second_identity=round(second[1][0], 4) if len(second) > 1 else None))
    if n % 100 == 0: print(n, aid[:8], r, round(ident, 3), flush=True)

df = pd.DataFrame(rows).sort_values("reads", ascending=False)
# confidence tiers: >=97% species-level, 94–97% genus-level plausible, <94% leave as "Prevotellaceae_unresolved"
df["call"] = np.select([df.identity >= 0.97, df.identity >= 0.94], [df.new_genus, df.new_genus + "_lowconf"], "Prevotellaceae_unresolved")
df.to_csv(OUT, index=False)
print(df.groupby("call").agg(n_asv=("reads", "size"), reads=("reads", "sum")).sort_values("reads", ascending=False))
print(df.groupby(["new_genus", "species"]).reads.sum().sort_values(ascending=False).head(20))
