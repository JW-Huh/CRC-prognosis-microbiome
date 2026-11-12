"""Rebuild the handover panels P0-P7 (SVG + PNG) from the tables in 4-output/.
Run from 5-Reanalysis-2026/:  python 2-script/make_panels2.py
"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pandas as pd, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator
from scipy.stats import spearmanr
from cox import coxph, logrank

SEG, LEY, PRV, GRY, GRN = "#2a78d6", "#eb6834", "#eda100", "#52514e", "#1baf7a"
plt.rcParams.update({"svg.fonttype": "none", "font.family": "DejaVu Sans", "font.size": 11,
                     "axes.spines.top": False, "axes.spines.right": False,
                     "axes.linewidth": 0.8, "figure.facecolor": "white", "axes.labelsize": 11})
OUT = "3-figure/"
def save(fig, n):
    fig.savefig(OUT + n + ".svg", bbox_inches="tight"); fig.savefig(OUT + n + ".png", dpi=200, bbox_inches="tight")
    plt.close(fig); print("ok", n)
def it(s): return r"$\it{" + s.replace(" ", r"\ ") + "}$"

def km(t, e):
    t = np.asarray(t, float); e = np.asarray(e, float); o = np.argsort(t); t, e = t[o], e[o]
    S, surv, times = 1.0, [1.0], [0.0]
    for ut in np.unique(t[e == 1]):
        n = (t >= ut).sum(); d = ((t == ut) & (e == 1)).sum(); S *= 1 - d / n; times.append(ut); surv.append(S)
    return np.array(times), np.array(surv)

d = pd.read_csv("4-output/diet2022_115_newtax.csv", index_col=0)
d["f46_hi"] = (d.Food46 >= d.Food46.median()).astype(int)
mb = pd.read_csv("4-output/merged_333_newtax.csv", index_col=0)

# P1 crosswalk
re_ = pd.read_csv("4-output/prevotella_reassignment.csv"); re_["c2"] = re_.call.str.replace("_lowconf", "")
sh = re_.groupby("c2").reads.sum().sort_values(ascending=False); sh = sh / sh.sum() * 100
lab = {"Segatella": "Segatella\n(" + it("S. copri") + ")", "Leyella": "Leyella\n(" + it("L. stercorea") + ")",
       "Prevotella": "Prevotella\ns. str. (oral)", "Hallella": "Hallella", "Hoylesella": "Hoylesella",
       "Prevotellaceae_unresolved": "Unresolved\n(<94% id)"}
fig, ax = plt.subplots(figsize=(5.0, 3.3))
ax.bar(range(len(sh)), sh.values, color=[{"Segatella": SEG, "Leyella": LEY, "Prevotella": PRV}.get(g, "#b6bec2") for g in sh.index], width=0.6)
for i, v in enumerate(sh.values): ax.text(i, v + 1.5, f"{v:.1f}%", ha="center", fontsize=10)
ax.set_xticks(range(len(sh))); ax.set_xticklabels([lab.get(g, g) for g in sh.index], fontsize=8.5)
ax.set_ylabel("% of 2022 'Prevotella' reads"); ax.set_ylim(0, 96)
save(fig, "P1_crosswalk_reads")

# P2 KM Food46
fig, ax = plt.subplots(figsize=(4.4, 3.5))
for g, c, l in [(1, GRN, "High intake (n=57, 4 events)"), (0, GRY, "Low intake (n=57, 18 events)")]:
    s = d[d.f46_hi == g]; t, su = km(s.PFS_time, s.PFS)
    ax.step(np.append(t, s.PFS_time.max()), np.append(su, su[-1]), where="post", color=c, lw=2, label=l)
    cen = s[s.PFS == 0]; i = np.searchsorted(t, cen.PFS_time, side="right") - 1
    ax.plot(cen.PFS_time, su[i], "|", color=c, ms=7, mew=1)
ax.set_ylim(0, 1.03); ax.set_xlabel("Days after surgery"); ax.set_ylabel("Progression-free survival")
ax.legend(frameon=False, fontsize=9, loc="lower left", bbox_to_anchor=(0.0, 0.03))
ax.text(0.03, 0.26, "log-rank p = 0.00055\nHR 0.18 (0.06–0.52), adj. stage", transform=ax.transAxes, fontsize=9, color=GRY)
save(fig, "P2_KM_Food46_PFS_114")

# P3 forest
cx = pd.read_csv("4-output/cox_prevotella_split_333.csv")
rows = [("'Prevotella' (2022, combined)", "Prevotella_2022", GRY), ("Segatella", "Segatella", SEG),
        ("Leyella", "Leyella", LEY), ("Prevotella s. str. (oral)", "Prevotella", PRV)]
fig, axes = plt.subplots(1, 2, figsize=(7.4, 2.7), sharey=True)
for ax, oc in zip(axes, ["PFS", "OS"]):
    for y, (l, tx, c) in zip(np.arange(len(rows))[::-1], rows):
        r = cx[(cx.taxon == tx) & (cx.outcome == oc)].iloc[0]
        ax.plot([r.lo, r.hi], [y, y], color=c, lw=1.8); ax.plot(r.HR_hi, y, "o", color=c, ms=7)
        ax.annotate(f"{r.HR_hi:.2f}", (r.HR_hi, y), textcoords="offset points", xytext=(0, 8), ha="center", fontsize=8.5, color=GRY)
    ax.axvline(1, color="#aab2b6", lw=0.8, ls="--"); ax.set_xscale("log"); ax.set_xlim(0.15, 3.2)
    ax.xaxis.set_major_locator(FixedLocator([0.25, 0.5, 1, 2])); ax.xaxis.set_minor_locator(NullLocator())
    ax.set_xticklabels(["0.25", "0.5", "1", "2"]); ax.set_title(oc, fontsize=11)
    ax.set_xlabel("Hazard ratio (high vs low CLR)")
axes[0].set_yticks(np.arange(len(rows))[::-1]); axes[0].set_yticklabels([r[0] for r in rows], fontsize=9.5)
save(fig, "P3_forest_taxa_333")

# P4 scatter
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1), sharex=True)
rk = d.Food46.rank()
for ax, tx, c in zip(axes, ["Leyella", "Segatella"], [LEY, SEG]):
    rho, p = spearmanr(d.Food46, d[tx])
    ax.scatter(rk, d[tx], s=20, color=c, alpha=.75, edgecolors="white", linewidths=.4)
    b = np.polyfit(rk, d[tx], 1); x = np.linspace(1, len(d), 50)
    ax.plot(x, np.polyval(b, x), color="#6b7478", lw=1, ls="--")
    ax.set_xlabel("Food46 intake rank (1–114)"); ax.set_title(it(tx), fontsize=12, color=c)
    ax.text(.03, .97, f"Spearman rho = {rho:+.2f}\np = {p:.3g}", transform=ax.transAxes, va="top", fontsize=9.5)
axes[0].set_ylabel("CLR abundance")
save(fig, "P4_scatter_Food46_vs_taxa")

# P5 lollipop
R = pd.read_csv("4-output/food46_vs_genera_clr_115.csv").sort_values("rho")
top = pd.concat([R.head(6), R.tail(6)])
fig, ax = plt.subplots(figsize=(4.9, 3.5))
cols = [LEY if g == "Leyella" else (SEG if g == "Segatella" else "#b6bec2") for g in top.genus]
ys = np.arange(len(top))
ax.hlines(ys, 0, top.rho, color=cols, lw=1.6)
for y, (r, c) in enumerate(zip(top.rho, cols)): ax.plot(r, y, "o", ms=6, color=c)
ax.set_yticks(ys); ax.set_yticklabels([it(g) if g in ("Leyella", "Segatella") else g for g in top.genus], fontsize=8.5)
ax.axvline(0, color="#2b3134", lw=0.8)
ax.set_xlabel("Spearman rho with Food46"); ax.set_title("70 genera screened (prevalence ≥20%)", fontsize=10)
save(fig, "P5_lollipop_genuswide")

# P6 PICRUSt2 — keep only pathways with meaningful output in at least one genus
prof = pd.read_csv("4-output/picrust_profile_Seg_vs_Ley.csv", index_col=0)
prof = prof[(prof.Seg_frac > 0.002) | (prof.Ley_frac > 0.002)].copy()
prof["diff"] = prof.Seg_frac - prof.Ley_frac
sel = pd.concat([prof.nsmallest(4, "diff"), prof.nlargest(4, "diff")]).sort_values("diff")
nm = {"HISDEG-PWY": "L-histidine\ndegradation", "PWY-7456": "β-mannan\ndegradation",
      "RHAMCAT-PWY": "L-rhamnose\ncatabolism", "PWY-6507": "4-deoxy-L-threo-\nhex-4-enopyranuronate deg.",
      "PWY-7242": "D-fructuronate\ndegradation", "TRPSYN-PWY": "L-tryptophan\nbiosynthesis"}
fig, ax = plt.subplots(figsize=(6.0, 3.2))
x = np.arange(len(sel)); w = .38
ax.bar(x - w/2, sel.Seg_frac * 100, w, color=SEG, label=it("Segatella"))
ax.bar(x + w/2, sel.Ley_frac * 100, w, color=LEY, label=it("Leyella"))
ax.set_xticks(x); ax.set_xticklabels([nm.get(i, i) for i in sel.index], fontsize=7.5, rotation=30, ha="right")
ax.set_ylabel("% of genus' predicted\npathway output"); ax.legend(frameon=False, fontsize=9.5)
save(fig, "P6_picrust_Seg_vs_Ley")

# P7 KM enterotype
fig, ax = plt.subplots(figsize=(4.4, 3.5))
for et, c, l in [("Segatella", SEG, "Segatella-type (n=158)"), ("Bacteroides", GRY, "Bacteroides-type (n=175)")]:
    s = mb[mb.ET_new == et]; t, su = km(s.OS_time, s.OS)
    ax.step(np.append(t, s.OS_time.max()), np.append(su, su[-1]), where="post", color=c, lw=2, label=l)
    cen = s[s.OS == 0]; i = np.searchsorted(t, cen.OS_time, side="right") - 1
    ax.plot(cen.OS_time, su[i], "|", color=c, ms=6, mew=.9)
ax.set_ylim(.5, 1.03); ax.set_xlabel("Days after surgery"); ax.set_ylabel("Overall survival")
ax.legend(frameon=False, fontsize=9, loc="lower left")
ax.text(.97, .96, "HR 0.50 (0.25–0.98), p = 0.045", transform=ax.transAxes, ha="right", va="top", fontsize=9, color=GRY)
save(fig, "P7_KM_enterotype_OS_333")

# P0 volcano  (uses the screening table written by the main pipeline, script 06)
dr = pd.read_csv("../3-Output/06_diet_vs_prognosis.csv")
dr["lh"] = np.log(dr.PFS_HR.astype(float)); dr["mp"] = -np.log10(dr.PFS_p_logrank.astype(float))
fig, ax = plt.subplots(figsize=(4.7, 3.6))
ax.axhline(2, ls=":", color="#aab2b6", lw=.9); ax.axvline(0, ls=":", color="#aab2b6", lw=.9)
ax.scatter(dr.lh, dr.mp, s=24, c=np.clip(dr.lh, -1.5, 1.5), cmap="coolwarm", edgecolors="#6b7478", linewidths=.3)
nm2 = {"Food46": "Food46 (stem namul)", "var18": "Retinol", "var83": "Vitamin D", "var63": "Leucine", "var65": "Methionine"}
for _, r in dr.nlargest(5, "mp").iterrows():
    ax.annotate(nm2.get(r["var"], r["var"]), (r.lh, r.mp), textcoords="offset points", xytext=(7, 2), fontsize=8.5)
ax.set_xlabel("log hazard ratio (PFS)"); ax.set_ylabel("−log10 log-rank p")
save(fig, "P0_volcano_2022_screen")
print("all panels rebuilt")
