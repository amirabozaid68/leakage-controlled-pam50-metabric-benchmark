"""Manuscript figures from ML_RESULTS_v1.2 (PAM50 benchmark, METABRIC + SCAN-B).
Usage (Colab or local):  python make_figures_v1.2.py <ML_RESULTS_v1.2 folder> <output folder>
Writes Fig2-Fig5 (main) and FigS1-FigS3 (supplementary) as PDF + 600-dpi PNG."""
import os, sys
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = sys.argv[1] if len(sys.argv) > 1 else "ML_RESULTS_v1.2"
FIG = sys.argv[2] if len(sys.argv) > 2 else "figures_v1.2"
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 7.5, "ytick.labelsize": 7.5,
    "legend.fontsize": 7.5, "axes.spines.top": False, "axes.spines.right": False,
    "axes.edgecolor": "#52514e", "axes.linewidth": 0.6, "xtick.color": "#52514e", "ytick.color": "#52514e",
    "axes.labelcolor": "#0b0b0b", "text.color": "#0b0b0b", "axes.grid": True, "grid.color": "#e6e5e0",
    "grid.linewidth": 0.5, "axes.axisbelow": True, "pdf.fonttype": 42, "ps.fonttype": 42})
INK2 = "#52514e"
W2 = 7.09                                              # double-column width (180 mm)

ARMS = {  # key: (label, colour)
    "pam50_centroid_v1": ("Centroid\n(v1.0)", "#b5b4ad"),
    "pam50_genefu":      ("genefu\nPAM50", "#52514e"),
    "lr_pam50genes":     ("LR\nPAM50\ngenes", "#2a78d6"),
    "lr_full3000":       ("LR\n3,000\ngenes", "#eb6834"),
    "compact_TOP80":     ("LR\nTOP80\npanel", "#1baf7a")}
COH = {"METABRIC": "METABRIC (microarray)", "SCANB": "SCAN-B (RNA-seq)"}
CLASSES = ["Basal", "Her2", "LumA", "LumB"]

ps = pd.read_csv(f"{RES}/per_seed/per_seed_metrics_all.tsv", sep="\t")
preds = pd.read_csv(f"{RES}/per_seed/predictions_all.tsv.gz", sep="\t")

def save(fig, name):
    for ext, kw in (("pdf", {}), ("png", {"dpi": 600})):
        fig.savefig(f"{FIG}/{name}.{ext}", bbox_inches="tight", facecolor="white", **kw)
    plt.close(fig); print("saved", name)

def panel_tag(ax, t):
    ax.text(-0.14, 1.06, t, transform=ax.transAxes, fontsize=11, fontweight="bold", va="bottom")

def dots(ax, df, metric, arms, ylim):
    rng = np.random.default_rng(1)
    for i, a in enumerate(arms):
        v = df.loc[df.arm == a, metric].to_numpy()
        lab, col = ARMS[a]
        ax.scatter(i + rng.uniform(-0.12, 0.12, len(v)), v, s=14, color=col, alpha=0.55, lw=0, zorder=3)
        m, sd = v.mean(), v.std(ddof=1)
        ax.errorbar(i + 0.28, m, yerr=sd, fmt="o", ms=5, color=col, mec="white", mew=1, elinewidth=1.4, capsize=0, zorder=4)
        ax.text(i + 0.28, m + sd + 0.006, f"{m:.3f}", ha="center", va="bottom", fontsize=6.8, color=INK2)
    ax.set_xticks(range(len(arms))); ax.set_xticklabels([ARMS[a][0] for a in arms])
    ax.set_xlim(-0.6, len(arms) - 0.2); ax.set_ylim(*ylim); ax.grid(axis="x", visible=False)

# ---------------- Fig 2: within-cohort performance -----------------------------------------
arms5 = list(ARMS)
fig, axs = plt.subplots(1, 2, figsize=(W2, 2.9), sharey=True)
for ax, c, t in zip(axs, COH, "ab"):
    dots(ax, ps[(ps.scope == "within") & (ps.test_cohort == c)], "macro_f1", arms5, (0.64, 1.0))
    ax.set_title(f"Within {COH[c]}", loc="left"); panel_tag(ax, t)
axs[0].set_ylabel("Macro-F1 on held-out test set")
fig.text(0.5, -0.04, "Dots: five 70/15/15 splits (seeds 123–127). Circle and bar: mean ± SD.",
         ha="center", fontsize=7, color=INK2)
fig.tight_layout(); save(fig, "Fig2_within_cohort_macroF1")

# ---------------- Fig 3: transfer across cohorts ---------------------------------------------
sup = ["lr_pam50genes", "lr_full3000", "compact_TOP80"]
fig, axs = plt.subplots(1, 2, figsize=(W2, 2.8), sharey=True)
for ax, tgt, t in zip(axs, COH, "ab"):
    src = [c for c in COH if c != tgt][0]
    gf = ps[(ps.arm == "pam50_genefu") & (ps.scope == "transfer") & (ps.test_cohort == tgt)].macro_f1.iloc[0]
    ax.axhline(gf, color=ARMS["pam50_genefu"][1], lw=1, ls=(0, (4, 3)), zorder=1)
    ax.text(-0.55, 0.755, f"dashed: genefu PAM50, no training ({gf:.3f})",
            ha="left", va="bottom", fontsize=6.6, color=INK2)
    for i, a in enumerate(sup):
        w = ps[(ps.arm == a) & (ps.scope == "within") & (ps.test_cohort == tgt)].macro_f1
        x = ps[(ps.arm == a) & (ps.scope == "transfer") & (ps.test_cohort == tgt)].macro_f1
        col = ARMS[a][1]
        ax.plot([i - 0.17, i + 0.17], [w.mean(), x.mean()], color=col, lw=1.2, zorder=2)
        ax.errorbar(i - 0.17, w.mean(), yerr=w.std(ddof=1), fmt="o", ms=5, mfc="white", mec=col, mew=1.3,
                    color=col, elinewidth=1.2, zorder=3)
        ax.errorbar(i + 0.17, x.mean(), yerr=x.std(ddof=1), fmt="o", ms=5, color=col, mec="white", mew=1,
                    elinewidth=1.2, zorder=3)
        ax.text(i + 0.24, x.mean(), f"{x.mean() - w.mean():+.3f}", va="center", fontsize=6.8, color=INK2)
    ax.set_xticks(range(len(sup))); ax.set_xticklabels([ARMS[a][0] for a in sup])
    ax.set_xlim(-0.6, len(sup) - 0.3); ax.set_ylim(0.75, 0.97); ax.grid(axis="x", visible=False)
    ax.set_title(f"Tested on {COH[tgt]}", loc="left"); panel_tag(ax, t)
    ax.text(0.0, -0.40, f"open: trained within {tgt.replace('SCANB', 'SCAN-B')}   filled: trained on "
            f"{src.replace('SCANB', 'SCAN-B')}", transform=ax.transAxes, fontsize=6.8, color=INK2)
axs[0].set_ylabel("Macro-F1")
fig.tight_layout(); save(fig, "Fig3_cross_cohort_transfer")

# ---------------- Fig 4: luminal boundary ----------------------------------------------------
def lum_rates(q):
    a, b = q[q.y_true == 2], q[q.y_true == 3]
    return (a.y_pred == 3).mean() * 100, (b.y_pred == 2).mean() * 100

rows = []
for tgt in COH:
    src = [c for c in COH if c != tgt][0]
    for a in ["pam50_genefu", "lr_pam50genes", "lr_full3000", "compact_TOP80"]:
        scopes = [(f"{tgt}_within", "within")] if a == "pam50_genefu" else \
                 [(f"{tgt}_within", "within"), (f"{src}_to_{tgt}", "transfer")]
        for sc, kind in scopes:
            q = preds[(preds.scope == sc) & (preds.arm == a)]
            if a == "pam50_genefu":
                q = preds[(preds.scope == f"{src}_to_{tgt}") & (preds.arm == a)]   # full target cohort, no training
            ab, ba = lum_rates(q)
            rows.append({"target": tgt, "arm": a, "kind": kind, "LumA_to_LumB": ab, "LumB_to_LumA": ba})
lum = pd.DataFrame(rows)
lum.to_csv(f"{FIG}/Fig4_source_luminal_rates.tsv", sep="\t", index=False)

fig, axs = plt.subplots(1, 2, figsize=(W2, 3.1), sharex=True)
order = [("pam50_genefu", "within"), ("lr_pam50genes", "within"), ("lr_pam50genes", "transfer"),
         ("lr_full3000", "within"), ("lr_full3000", "transfer"), ("compact_TOP80", "within"), ("compact_TOP80", "transfer")]
for ax, tgt, t in zip(axs, COH, "ab"):
    src = [c for c in COH if c != tgt][0].replace("SCANB", "SCAN-B")
    d = lum[lum.target == tgt].set_index(["arm", "kind"])
    ys = np.arange(len(order))[::-1]
    for y, key in zip(ys, order):
        r = d.loc[key]; col = ARMS[key[0]][1]
        hatch = "////" if key[1] == "transfer" else None
        ax.barh(y + 0.19, -r.LumA_to_LumB, 0.36, color=col if hatch is None else "white", edgecolor=col, hatch=hatch, lw=0.8)
        ax.barh(y - 0.19, r.LumB_to_LumA, 0.36, color=col if hatch is None else "white", edgecolor=col, hatch=hatch, lw=0.8)
        ax.text(-r.LumA_to_LumB - 0.8, y + 0.19, f"{r.LumA_to_LumB:.1f}", ha="right", va="center", fontsize=6.3, color=INK2)
        ax.text(r.LumB_to_LumA + 0.8, y - 0.19, f"{r.LumB_to_LumA:.1f}", ha="left", va="center", fontsize=6.3, color=INK2)
    names = {"pam50_genefu": "genefu PAM50 (no training)", "lr_pam50genes": "LR, PAM50 genes",
             "lr_full3000": "LR, 3,000 genes", "compact_TOP80": "LR, TOP80 panel"}
    labels = [names[a] + ("" if a == "pam50_genefu" else (", other cohort" if k == "transfer" else ", same cohort"))
              for a, k in order]
    ax.set_yticks(ys); ax.set_yticklabels(labels if t == "a" else [])
    ax.axvline(0, color=INK2, lw=0.6); ax.grid(axis="y", visible=False)
    ax.set_xlim(-42, 42); ax.set_xticks([-30, -20, -10, 0, 10, 20, 30]); ax.set_xticklabels(["30", "20", "10", "0", "10", "20", "30"])
    ax.set_title(f"Tested on {COH[tgt]}", loc="left", pad=14); panel_tag(ax, t)
    ax.text(0.02, 1.0, "LumA called LumB", transform=ax.transAxes, fontsize=6.8, color=INK2, va="bottom")
    ax.text(0.98, 1.0, "LumB called LumA", transform=ax.transAxes, fontsize=6.8, color=INK2, va="bottom", ha="right")
from matplotlib.patches import Patch
fig.legend(handles=[Patch(facecolor=INK2, edgecolor=INK2, label="model trained in the tested cohort (held-out test sets)"),
                    Patch(facecolor="white", edgecolor=INK2, hatch="////", label="model trained in the other cohort (every 4-class tumour)")],
           loc="lower center", ncol=2, frameon=False, bbox_to_anchor=(0.55, -0.10))
fig.text(0.62, 0.0, "Tumours assigned to the other luminal subtype (%)", ha="center", fontsize=8)
fig.tight_layout(w_pad=1.0); save(fig, "Fig4_luminal_boundary")

# ---------------- Fig 5: panel stability and panel size ---------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(W2, 2.5), gridspec_kw={"width_ratios": [1.1, 0.9, 1.3]})
ax = axs[0]
for j, (c, col) in enumerate(zip(COH, ["#2a78d6", "#eb6834"])):
    f = pd.read_csv(f"{RES}/stability/gene_inclusion_frequency_{c}.tsv", sep="\t")
    n = f.seed_count.value_counts().reindex(range(1, 6), fill_value=0)
    ax.bar(n.index + (j - 0.5) * 0.38, n.values, 0.36, color=col, label=c.replace("SCANB", "SCAN-B"))
    for k, v in n.items(): ax.text(k + (j - 0.5) * 0.38, v + 2, str(v), ha="center", fontsize=6.3, color=INK2)
ax.set_xlabel("Splits in which a gene entered TOP80"); ax.set_ylabel("Genes"); ax.set_xticks(range(1, 6))
ax.legend(frameon=False, loc="upper right"); ax.grid(axis="x", visible=False); panel_tag(ax, "a")

ax = axs[1]
st = pd.read_csv(f"{RES}/stability/stability_summary.tsv", sep="\t").set_index("cohort")
a, b = st.loc["METABRIC", "stable_size"], st.loc["SCANB", "stable_size"]
ov = st.loc["METABRIC_and_SCANB_stable_overlap", "stable_size"]
ov_p = st.loc["METABRIC_and_SCANB_stable_overlap", "stable_in_PAM50"]
vals = [(f"METABRIC only", a - ov, "#2a78d6"), ("Both cohorts", ov, INK2), ("SCAN-B only", b - ov, "#eb6834")]
for i, (lab, v, col) in enumerate(vals):
    ax.bar(i, v, 0.6, color=col); ax.text(i, v + 1, str(int(v)), ha="center", fontsize=7, color=INK2)
ax.set_xticks(range(3)); ax.set_xticklabels([v[0].replace(" ", "\n", 1) for v in vals], fontsize=6.6)
ax.set_ylabel("Stable genes (≥3 of 5 splits)"); ax.grid(axis="x", visible=False); panel_tag(ax, "b")
ax.text(1, ov + 6, f"{int(ov_p)} in PAM50", ha="center", fontsize=6.5, color=INK2)

ax = axs[2]
ks = [20, 40, 50, 80, 120]
for c, col, mk in zip(COH, ["#2a78d6", "#eb6834"], ["o", "s"]):
    w = ps[(ps.scope == "within") & (ps.test_cohort == c)]
    m = [w[w.arm == f"compact_TOP{k}"].macro_f1.mean() for k in ks]
    s = [w[w.arm == f"compact_TOP{k}"].macro_f1.std(ddof=1) for k in ks]
    ax.errorbar(range(len(ks)), m, yerr=s, color=col, marker=mk, ms=4, lw=1.4, elinewidth=0.9, capsize=0,
                label=f"{c.replace('SCANB', 'SCAN-B')}, compact panel (LR)")
    ref = w[w.arm == "lr_pam50genes"].macro_f1.mean()
    ax.axhline(ref, color=col, lw=0.9, ls=(0, (4, 3)), label=f"{c.replace('SCANB', 'SCAN-B')}, LR on 50 PAM50 genes")
ax.set_xlabel("Compact panel size (genes)"); ax.set_ylabel("Macro-F1 (within cohort)")
ax.set_xticks(range(len(ks))); ax.set_xticklabels(ks); ax.set_xlim(-0.4, len(ks) - 0.6)
ax.set_ylim(0.64, 0.96); ax.legend(frameon=False, loc="lower right", ncol=1, fontsize=5.8, handlelength=1.8)
panel_tag(ax, "c")
fig.tight_layout(w_pad=1.2); save(fig, "Fig5_panel_stability_and_size")

# ---------------- Fig S1: all four metrics, all arms -----------------------------------------
SHORT = ["Centroid v1.0", "genefu", "LR PAM50 genes", "LR 3,000 genes", "LR TOP80"]
mets = [("macro_f1", "Macro-F1"), ("bacc", "Balanced accuracy"), ("macro_roc_auc", "Macro ROC-AUC"), ("macro_pr_auc", "Macro PR-AUC")]
fig, axs = plt.subplots(2, 4, figsize=(W2, 4.8))
for r, c in enumerate(COH):
    w = ps[(ps.scope == "within") & (ps.test_cohort == c)]
    for k, (m, lab) in enumerate(mets):
        ax = axs[r, k]
        lo = max(0.6, w[m].min() - 0.03)
        for i, a in enumerate(arms5):
            v = w.loc[w.arm == a, m]
            ax.errorbar(i, v.mean(), yerr=v.std(ddof=1), fmt="o", ms=4, color=ARMS[a][1], elinewidth=1, capsize=0)
        ax.set_xticks(range(5)); ax.set_xticklabels(SHORT, fontsize=6.3, rotation=45, ha="right", rotation_mode="anchor")
        ax.set_ylim(lo, 1.0); ax.grid(axis="x", visible=False)
        if r == 0: ax.set_title(lab)
        if k == 0: ax.set_ylabel(c.replace("SCANB", "SCAN-B"))
fig.tight_layout(); save(fig, "FigS1_all_metrics_within")

# ---------------- Fig S2: per-class recall heatmap -------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(W2, 2.6))
for ax, c, t in zip(axs, COH, "ab"):
    w = ps[(ps.scope == "within") & (ps.test_cohort == c)]
    M = np.array([[w[w.arm == a][f"recall_{k}"].mean() for k in CLASSES] for a in arms5])
    im = ax.imshow(M, cmap="Blues", vmin=0.5, vmax=1.0, aspect="auto")
    for i in range(M.shape[0]):
        for j in range(4):
            ax.text(j, i, f"{M[i, j]:.2f}", ha="center", va="center", fontsize=7, color="white" if M[i, j] > 0.8 else "#0b0b0b")
    ax.set_xticks(range(4)); ax.set_xticklabels(CLASSES); ax.set_yticks(range(5))
    ax.set_yticklabels([ARMS[a][0].replace("\n", " ") for a in arms5] if t == "a" else []); ax.grid(False)
    ax.set_title(f"Recall within {c.replace('SCANB', 'SCAN-B')}", loc="left"); panel_tag(ax, t)
fig.colorbar(im, ax=axs, shrink=0.8, label="Mean recall over 5 splits")
save(fig, "FigS2_per_class_recall")

# ---------------- Fig S3: genefu against published labels ------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(W2, 2.9))
for ax, c, t in zip(axs, COH, "ab"):
    ct = pd.read_csv(f"{RES}/genefu/genefu_vs_published_crosstab_{c}.tsv", sep="\t", index_col=0)
    R = ct.div(ct.sum(axis=1), axis=0)
    ax.imshow(R.values, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            ax.text(j, i, f"{ct.values[i, j]}", ha="center", va="center", fontsize=6.8,
                    color="white" if R.values[i, j] > 0.6 else "#0b0b0b")
    ax.set_xticks(range(R.shape[1])); ax.set_xticklabels(R.columns); ax.set_yticks(range(R.shape[0]))
    ax.set_yticklabels(R.index); ax.set_xlabel("genefu call"); ax.set_ylabel("Published PAM50 label"); ax.grid(False)
    ax.set_title(COH[c], loc="left"); panel_tag(ax, t)
fig.tight_layout(); save(fig, "FigS3_genefu_vs_published_labels")

# ---------------- Fig 1: study design ---------------------------------------------------------
from matplotlib.patches import FancyBboxPatch
fig, ax = plt.subplots(figsize=(W2, 5.8)); ax.set_xlim(0, 100); ax.set_ylim(0, 82); ax.axis("off")
def box(x, y, w, h, title, body, edge=INK2):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.2,rounding_size=1.2", fc="#f7f7f5", ec=edge, lw=0.9))
    ax.text(x + 1.5, y + h - 1.3, title, fontsize=7.8, fontweight="bold", va="top")
    ax.text(x + 1.5, y + h - 4.8, body, fontsize=6.4, va="top", color="#2b2b29", linespacing=1.4)
def arrow(x0, y0, x1, y1):
    ax.annotate("", (x1, y1), (x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=0.9, shrinkA=0, shrinkB=0))
box(2, 67, 45, 14, "METABRIC (Illumina microarray)",
    "1,981 profiles; 1,974 with a PAM50 label\n4-class tumours: n = 1,608\nBasal 209 · Her2 224 · LumA 700 · LumB 475", edge="#2a78d6")
box(53, 67, 45, 14, "SCAN-B, GSE96058 (RNA-seq)",
    "3,273 primary tumours (replicates removed)\n4-class tumours: n = 3,052\nBasal 339 · Her2 327 · LumA 1,657 · LumB 729", edge="#eb6834")
box(28, 53.5, 44, 10, "Shared gene space", "genes measured in both cohorts\ncurrent PAM50 gene symbols")
arrow(24.5, 67, 36, 63.4); arrow(75.5, 67, 64, 63.4)
box(2, 17, 45, 34, "A  Within-cohort benchmark",
    "5 seeds, stratified 70/15/15 splits\n"
    "preprocessing fitted on the training set only\n"
    "C tuned on validation (macro-F1),\nrefit on train + validation, one test pass\n\n"
    "Arms on identical splits\n"
    "• genefu PAM50 (test set alone)\n"
    "• LR on the 50 PAM50 genes\n"
    "• LR on 3,000 train-filtered genes\n"
    "• L1-selected panel + LR (20–120 genes)\n"
    "• v1.0 centroid re-implementation")
box(53, 17, 45, 34, "B  Cross-cohort transfer",
    "both directions:\nMETABRIC → SCAN-B and SCAN-B → METABRIC\n\n"
    "every model refit in A predicts every\n4-class tumour of the other cohort\n\n"
    "target cohort z-scored on itself,\nlabels never used\n\n"
    "reference: genefu PAM50 on the whole\ntarget cohort, no training")
arrow(42, 53.5, 26, 51.4); arrow(58, 53.5, 74, 51.4)
box(2, 2.5, 96, 11, "Outputs",
    "macro-F1, balanced accuracy, macro ROC-AUC and PR-AUC · paired bootstrap CI and McNemar test per split\n"
    "LumA ↔ LumB reassignment · genefu agreement with published labels · panel stability and cross-cohort overlap")
arrow(24.5, 17, 24.5, 13.9); arrow(75.5, 17, 75.5, 13.9)
save(fig, "Fig1_study_design")
