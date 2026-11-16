#!/usr/bin/env python3
"""
Micro-derived chromosome composition analysis.

Panels:
  A  Per-species scatter: log10(chr size Mb) vs GC content  — all species in a grid
  B  Per-species scatter: log10(chr size Mb) vs gene density
  C  Within-lacertid strip plot: macro vs micro-derived vs true-micro (GC + gene density)
  D  Outgroup-only: GC of micro-derived chromosomes vs divergence time from Lacertidae crown
  E  Same as D but gene density

Key distinctions:
  - "true micro"    : Chr18 (#157AB7) — never macrochromosomized
  - "micro-derived" : Chr2,8,13,14,15,16,ChrZ — ancestral micros now embedded in macros
  - "macro-derived" : grey (#C0C0C0) + Chr17 (#B2D78A)

Inputs:
  TSV       tab-separated table: species, chr, label, psa_color, size_mb, gc, gene_density
  TREE_FILE time-calibrated squamate species tree in Newick format (ultrametric)

Usage:
  python micro_decay_analysis.py

Requires: pandas, numpy, matplotlib, scipy, ete3
"""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats
from ete3 import Tree

# ── paths ─────────────────────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TSV        = os.path.join(SCRIPT_DIR, "outgroups_composition.tsv")
TREE_FILE  = os.path.join(SCRIPT_DIR, "squamate_species_tree.tree")
OUTDIR     = SCRIPT_DIR

# ── species metadata ──────────────────────────────────────────────────────────
ALL_SPECIES = [
    ("Gallus_gallus",               "G. gallus",        "outgroup"),
    ("Dibamus_smithi",              "D. smithi",        "outgroup"),
    ("Iguana_iguana",               "I. iguana",        "outgroup"),
    ("Varanus_acanthurus",          "V. acanthurus",    "outgroup"),
    ("Vipera_latastei",             "V. latastei",      "outgroup"),
    ("Aspidocelis_tigris",          "A. tigris",        "outgroup"),
    ("Rhineura_floridana",          "R. floridana",     "outgroup"),
    ("Psammodromus_algirus",        "P. algirus",       "lacertid"),
    ("Acanthodactylus_erythrurus",  "A. erythrurus",    "lacertid"),
    ("Eremias_argus",               "E. argus",         "lacertid"),
    ("Timon_lepidus",               "T. lepidus",       "lacertid"),
    ("Lacerta_agilis",              "L. agilis",        "lacertid"),
    ("Iberolacerta_aurelioi",       "I. aurelioi",      "lacertid"),
    ("Zootoca_vivipara",            "Z. vivipara",      "lacertid"),
    ("Scelarcis_perspicillata",     "S. perspicillata", "lacertid"),
    ("Podarcis_muralis",            "P. muralis",       "lacertid"),
]
OUTGROUP_SPS = [s for s, _, k in ALL_SPECIES if k == "outgroup"]
LACERTID_SPS = [s for s, _, k in ALL_SPECIES if k == "lacertid"]
SP_LABEL     = {s: l for s, l, _ in ALL_SPECIES}

# ── color / category scheme ───────────────────────────────────────────────────
COLOR_META = {
    "#C0C0C0": ("macro",       "Macro",    "#C0C0C0"),
    "#B2D78A": ("macro",       "Macro",    "#C0C0C0"),
    "#E96456": ("micro_mac",   "Chr2",     "#E96456"),
    "#FDBC6B": ("micro_mac",   "Chr8",     "#FDBC6B"),
    "#42956B": ("micro_mac",   "Chr13",    "#42956B"),
    "#BB8F6E": ("micro_mac",   "Chr14",    "#BB8F6E"),
    "#F6979A": ("micro_mac",   "Chr15",    "#F6979A"),
    "#8ED0C4": ("micro_mac",   "Chr16",    "#8ED0C4"),
    "#E36EAA": ("micro_mac",   "ChrZ",     "#E36EAA"),
    "#157AB7": ("true_micro",  "Chr18\n(true micro)", "#157AB7"),
}
PLOT_COLOR = {k: v[2] for k, v in COLOR_META.items()}

# ── load data ─────────────────────────────────────────────────────────────────
df = pd.read_csv(TSV, sep="\t")
df["category"] = df["psa_color"].map(lambda c: COLOR_META.get(c, ("macro","","#C0C0C0"))[0])
df["plot_col"]  = df["psa_color"].map(PLOT_COLOR)
df["kind"]      = df["species"].map(lambda s: "lacertid" if s in LACERTID_SPS else "outgroup")

# ── divergence times from Lacertidae crown ───────────────────────────────────
print("Computing divergence times from Lacertidae crown...")
t = Tree(TREE_FILE)
t.prune([s for s, _, _ in ALL_SPECIES], preserve_branch_length=True)

lacertid_nodes = [t.search_nodes(name=s)[0] for s in LACERTID_SPS]
lac_crown      = t.get_common_ancestor(lacertid_nodes)

def node_age(node):
    """Age of node = distance from node to any descendant leaf (ultrametric tree)."""
    leaf = next(node.iter_leaves())
    d, n = 0.0, leaf
    while n != node:
        d += n.dist
        n = n.up
    return d

lac_crown_age = node_age(lac_crown)
print(f"  Lacertidae crown age: {lac_crown_age:.1f} Mya")

div_times = {}
for sp in OUTGROUP_SPS:
    sp_node  = t.search_nodes(name=sp)[0]
    mrca     = t.get_common_ancestor([sp_node, lac_crown])
    div_times[sp] = node_age(mrca)
    print(f"  {sp}: {div_times[sp]:.1f} Mya from Lacertidae crown")

df["div_time"] = df["species"].map(div_times)

# ── helpers ───────────────────────────────────────────────────────────────────
def despine(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"

# ── FIGURE ────────────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(22, 28))

outer = gridspec.GridSpec(4, 1, figure=fig,
                          height_ratios=[3, 3, 2, 2],
                          hspace=0.38)

# ── A & B: per-species scatter grids ─────────────────────────────────────────
for row_idx, (y_key, y_lab, panel_letter) in enumerate([
        ("gc",           "GC content",             "A"),
        ("gene_density", "Gene density (genes/Mb)", "B"),
]):
    inner = gridspec.GridSpecFromSubplotSpec(
        4, 4, subplot_spec=outer[row_idx], hspace=0.55, wspace=0.35)

    for i, (sp, sp_lab, kind) in enumerate(ALL_SPECIES):
        ax = fig.add_subplot(inner[i // 4, i % 4])
        sub = df[df["species"] == sp]

        x = np.log10(sub["size_mb"].values)
        y = sub[y_key].values
        c = sub["plot_col"].values

        ax.scatter(x, y, c=c, edgecolors="white", linewidths=0.3, s=30, zorder=3)
        for xi, yi, lb in zip(x, y, sub["label"].values):
            ax.annotate(lb, (xi, yi), textcoords="offset points",
                        xytext=(2, 2), fontsize=4.5, color="#444")

        bg = "#eef4ee" if kind == "lacertid" else "#f8f8f8"
        ax.set_facecolor(bg)
        ax.set_title(rf"$\it{{{sp_lab}}}$", fontsize=7, pad=2)
        ax.set_xlabel("log10 size (Mb)", fontsize=6)
        ax.set_ylabel(y_lab, fontsize=6)
        ax.tick_params(labelsize=5.5)
        ax.grid(color="lightgrey", lw=0.4, alpha=0.6)
        despine(ax)

    fig.text(0.01, [0.985, 0.665][row_idx], panel_letter,
             fontsize=14, fontweight="bold", va="top")

# ── C: within-lacertid strip (GC and gene density) ───────────────────────────
lac_df = df[df["kind"] == "lacertid"].copy()

inner_c = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[2], wspace=0.30)

CAT_ORDER  = ["macro", "micro_mac", "true_micro"]
CAT_COLORS = {"macro": "#C0C0C0", "micro_mac": "#888888", "true_micro": "#157AB7"}
CAT_LABELS = {"macro": "Macro-derived", "micro_mac": "Micro-derived\n(in macro)", "true_micro": "True micro\n(Chr18)"}

rng = np.random.default_rng(0)

for ci, (y_key, y_lab) in enumerate([("gc", "GC content"),
                                      ("gene_density", "Gene density (genes/Mb)")]):
    ax = fig.add_subplot(inner_c[0, ci])

    positions = {c: i for i, c in enumerate(CAT_ORDER)}
    for cat in CAT_ORDER:
        sub = lac_df[lac_df["category"] == cat][y_key].dropna().values
        if len(sub) == 0:
            continue
        xi = positions[cat]
        jx = xi + rng.uniform(-0.18, 0.18, len(sub))
        dot_cols = lac_df[lac_df["category"] == cat]["plot_col"].values
        ax.scatter(jx, sub, c=dot_cols, edgecolors="white",
                   linewidths=0.3, s=28, alpha=0.85, zorder=3)
        ax.plot([xi - 0.22, xi + 0.22], [np.median(sub)] * 2,
                color="#222", lw=2, zorder=4)

    v_mac = lac_df[lac_df["category"] == "macro"][y_key].dropna().values
    v_mic = lac_df[lac_df["category"] == "micro_mac"][y_key].dropna().values
    v_tru = lac_df[lac_df["category"] == "true_micro"][y_key].dropna().values

    _, p1 = stats.mannwhitneyu(v_mac, v_mic, alternative="two-sided")
    _, p2 = stats.mannwhitneyu(v_mic, v_tru, alternative="two-sided")

    ymax = lac_df[y_key].max() * 1.04
    ax.annotate("", xy=(1, ymax), xytext=(0, ymax),
                arrowprops=dict(arrowstyle="-", color="#555"))
    ax.text(0.5, ymax * 1.005, sig_label(p1), ha="center", fontsize=8, color="#555")
    ax.annotate("", xy=(2, ymax * 1.07), xytext=(1, ymax * 1.07),
                arrowprops=dict(arrowstyle="-", color="#555"))
    ax.text(1.5, ymax * 1.075, sig_label(p2), ha="center", fontsize=8, color="#555")

    ax.set_xticks(range(len(CAT_ORDER)))
    ax.set_xticklabels([CAT_LABELS[c] for c in CAT_ORDER], fontsize=8)
    ax.set_ylabel(y_lab, fontsize=9)
    ax.tick_params(labelsize=7.5)
    ax.grid(axis="y", color="lightgrey", lw=0.5, alpha=0.6)
    ax.set_facecolor("#fafafa")
    despine(ax)

fig.text(0.01, 0.355, "C", fontsize=14, fontweight="bold", va="top")

# ── D & E: GC / gene density vs divergence time ───────────────────────────────
out_df    = df[df["kind"] == "outgroup"].copy()
out_micro = out_df[out_df["category"].isin(["micro_mac", "true_micro"])].copy()

lac_medians = (lac_df.groupby("category")[["gc", "gene_density"]].median())

inner_d = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[3], wspace=0.30)

for ci, (y_key, y_lab) in enumerate([("gc", "GC content"),
                                      ("gene_density", "Gene density (genes/Mb)")]):
    ax = fig.add_subplot(inner_d[0, ci])

    x = out_micro["div_time"].values
    y = out_micro[y_key].values
    c = out_micro["plot_col"].values

    ax.scatter(x, y, c=c, edgecolors="white", linewidths=0.4, s=45, zorder=3, alpha=0.85)

    for sp in OUTGROUP_SPS:
        sp_x = div_times.get(sp)
        if sp_x:
            ax.axvline(sp_x, color="lightgrey", lw=0.5, zorder=0)
            ax.text(sp_x, y.min() * 0.997,
                    SP_LABEL[sp].split(".")[1].strip(),
                    ha="center", va="top", fontsize=5.5, rotation=45, color="#555")

    for cat, col in [("micro_mac", "#888888"), ("true_micro", "#157AB7")]:
        if cat in lac_medians.index:
            med = lac_medians.loc[cat, y_key]
            ax.axhline(med, color=col, lw=1.2, ls="--", alpha=0.7,
                       label=f"Lacertid {CAT_LABELS[cat]} median")

    sub_reg = out_micro[out_micro["category"] == "micro_mac"]
    if len(sub_reg) >= 4:
        xr = sub_reg["div_time"].values
        yr = sub_reg[y_key].values
        slope, intercept, r, p, _ = stats.linregress(xr, yr)
        xl = np.linspace(xr.min(), xr.max(), 100)
        ax.plot(xl, slope * xl + intercept, color="#333", lw=1.5, alpha=0.7)
        ax.text(0.97, 0.97, f"r={r:.2f}, {sig_label(p)}",
                transform=ax.transAxes, ha="right", va="top", fontsize=7, color="#333")

    ax.axvspan(0, lac_crown_age, alpha=0.06, color="green",
               label=f"Lacertidae crown ({lac_crown_age:.0f} Mya)")
    ax.text(lac_crown_age / 2, y.max(),
            "Lacertidae", ha="center", fontsize=6.5, color="green", alpha=0.7)

    ax.set_xlabel("Divergence from Lacertidae crown (Mya)", fontsize=9)
    ax.set_ylabel(y_lab, fontsize=9)
    ax.tick_params(labelsize=7.5)
    ax.grid(color="lightgrey", lw=0.4, alpha=0.6)
    ax.set_facecolor("#fafafa")
    ax.legend(fontsize=6, loc="lower right", framealpha=0.8)
    despine(ax)

fig.text(0.01, 0.185, "D", fontsize=14, fontweight="bold", va="top")

# ── global legend ─────────────────────────────────────────────────────────────
CHROM_LEGEND = [
    ("#E96456", "Chr2"),  ("#FDBC6B", "Chr8"),  ("#42956B", "Chr13"),
    ("#BB8F6E", "Chr14"), ("#F6979A", "Chr15"),  ("#8ED0C4", "Chr16"),
    ("#E36EAA", "ChrZ"),  ("#157AB7", "Chr18 (true micro)"),
    ("#C0C0C0", "Macro (incl. Chr17)"),
]
handles = [plt.scatter([], [], c=col, s=30, label=lab, edgecolors="white")
           for col, lab in CHROM_LEGEND]
fig.legend(handles=handles, title="P. algirus chr",
           loc="lower center", ncol=5, fontsize=7, title_fontsize=7.5,
           bbox_to_anchor=(0.5, -0.01), framealpha=0.9)

fig.suptitle("Micro-derived chromosome composition across squamates",
             fontsize=13, y=1.002, fontweight="bold")

out_pdf = os.path.join(OUTDIR, "micro_decay_analysis.pdf")
fig.savefig(out_pdf, format="pdf", bbox_inches="tight")
plt.close(fig)
print(f"\nSaved: {out_pdf}")
