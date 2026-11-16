#!/usr/bin/env python3
"""
Enrichment-based GC content and gene density comparison panels.

Two variables are kept separate:
  - psa_color: homology to P. algirus chromosomes (synteny crosswalk)
  - is_macro : karyotype-based macro/micro status from chromosome_classification.txt

Outgroup reference (C–D): ONLY remaining macrochromosomes (macros that are not
  Chr2 or Chr13 homologues).  Outgroup micros are EXCLUDED from the reference
  even when their lacertid homologue is now a macrochromosome (e.g. Chr8 homologs
  are micro in every outgroup and so belong to "retained micro", not the reference).

Lacertid reference (G–H): remaining macrochromosomes (psa_color in grey/Chr17 set),
  excluding Chr2, Chr13, the five established micro-derived chromosomes and Chr18.

Iberolacerta aurelioi: present in the TSV (13 chromosomes; 2n = 24 karyotype,
  substantially more fused than the conserved 2n = 36–38 studied here) and therefore
  excluded from the lacertid panels (not part of the 8-species analysis set).

Groups: assigned by cross-referencing psa_color (homology) and macro/micro status
  (karyotype).  No conflicts detected: all Chr2 and Chr13 homologues are
  macrochromosomes in every outgroup species.

Inputs:
  TSV      tab-separated table with columns: species, chr, label, psa_color,
           size_mb, gc, gene_density
  CLF_TXT  plain-text chromosome classification file (see below)

Usage:
  python plot_gc_gene_density.py

Outputs:
  plot_gc_gene_density.pdf / .png
  chr_classification.tsv
  species_summary.tsv
"""

import os, re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.lines as mlines
from matplotlib.ticker import FormatStrFormatter

matplotlib.rcParams["font.family"] = ["Optima", "URW Classico", "Palatino", "DejaVu Sans"]

# ── input / output paths ───────────────────────────────────────────────────────
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TSV     = os.path.join(SCRIPT_DIR, "outgroups_composition.tsv")
CLF_TXT = os.path.join(SCRIPT_DIR, "chromosome_classification.txt")
OUTDIR  = SCRIPT_DIR

# ── species sets ───────────────────────────────────────────────────────────────
LACERTID_SP = {
    "Psammodromus_algirus", "Acanthodactylus_erythrurus", "Eremias_argus",
    "Timon_lepidus", "Lacerta_agilis", "Zootoca_vivipara",
    "Scelarcis_perspicillata", "Podarcis_muralis",
}
ALL_OUTGROUP_SP = [
    "Gallus_gallus", "Dibamus_smithi",
    "Vipera_latastei", "Iguana_iguana", "Varanus_acanthurus",
    "Rhineura_floridana", "Aspidocelis_tigris",
]
ENRICH_OUT_SP = [sp for sp in ALL_OUTGROUP_SP if sp != "Gallus_gallus"]

# ── species colours ────────────────────────────────────────────────────────────
LAC_COL        = "#44B2A8"
GREY_DARK      = "#555555"
GREY_MID       = "#ABABAB"
TOXICOFERA_COL = "#AA66A5"
RHINEURA_COL   = "#7492CB"
ASPIDO_COL     = "#EDE756"

SP_COLOR = {
    "Gallus_gallus":      GREY_DARK,
    "Dibamus_smithi":     GREY_MID,
    "Vipera_latastei":    TOXICOFERA_COL,
    "Iguana_iguana":      TOXICOFERA_COL,
    "Varanus_acanthurus": TOXICOFERA_COL,
    "Rhineura_floridana": RHINEURA_COL,
    "Aspidocelis_tigris": ASPIDO_COL,
    **{sp: LAC_COL for sp in LACERTID_SP},
}

LAC_SP_COL = {
    "Psammodromus_algirus":       "#1B7837",
    "Acanthodactylus_erythrurus": "#762A83",
    "Eremias_argus":              "#4393C3",
    "Timon_lepidus":              "#D6604D",
    "Lacerta_agilis":             "#F4A582",
    "Zootoca_vivipara":           "#92C5DE",
    "Scelarcis_perspicillata":    "#878787",
    "Podarcis_muralis":           "#E08214",
}

# ── chromosome-colour constants ────────────────────────────────────────────────
MACRO_PSA      = "#C0C0C0"
CHR2_PSA       = "#E96456"
CHR13_PSA      = "#42956B"
CHR17_PSA      = "#B2D78A"
CHR18_PSA      = "#157AB7"
MICRO_DER_COLS = frozenset({"#FDBC6B", "#BB8F6E", "#F6979A", "#8ED0C4", "#E36EAA"})
MACRO_REM_COLS = frozenset({MACRO_PSA, CHR17_PSA})

MACRO_FILL      = "#D0D0D0"
MICRO_DER_FILL  = "#F68C1F"
TRUE_MICRO_FILL = "#1477A8"

# ── 1. Parse macro/micro classification ───────────────────────────────────────
def _parse_classification(txt_path):
    result = {}; current = None; reading_micro = False
    with open(txt_path) as f:
        for raw in f:
            line = raw.rstrip("\n"); s = line.strip()
            if not s or s.startswith("="): reading_micro = False; continue
            indent = len(line) - len(line.lstrip(" "))
            if indent == 0:
                current = s; result[current] = {"macro": set(), "micro": set()}; reading_micro = False
            elif indent == 2 and "Macrochromosomes" in s:
                after = s.split(":", 1)[1] if ":" in s else ""
                result[current]["macro"] = {re.sub(r"^chr", "", x.strip())
                                             for x in after.split(",") if x.strip()}
                reading_micro = False
            elif indent == 2 and "Microchromosomes" in s:
                reading_micro = True
            elif indent == 4 and reading_micro:
                result[current]["micro"].update(x.strip() for x in s.split(",") if x.strip())
    return result

clf = _parse_classification(CLF_TXT)

def _get_macro_set(species_underscored):
    sp = species_underscored.replace("_", " ")
    return clf.get(sp, {}).get("macro", set())

def _get_micro_set(species_underscored):
    sp = species_underscored.replace("_", " ")
    return clf.get(sp, {}).get("micro", set())

# ── 2. Load & filter data ─────────────────────────────────────────────────────
df_raw = pd.read_csv(TSV, sep="\t")
all_sp  = set(SP_COLOR.keys())
df      = df_raw[df_raw["species"].isin(all_sp)].copy()
df["log_size"] = np.log10(df["size_mb"])

df_lac = df[df["species"].isin(LACERTID_SP)].copy()
df_out = df[df["species"].isin(ENRICH_OUT_SP)].copy()

# ── 3. Assign is_macro to outgroup chromosomes ─────────────────────────────────
def _is_macro_row(row):
    lbl = str(row["label"])
    mac = _get_macro_set(row["species"])
    mic = _get_micro_set(row["species"])
    if lbl in mac: return True
    if lbl in mic: return False
    return None

df_out["is_macro"] = df_out.apply(_is_macro_row, axis=1)

unc = df_out[df_out["is_macro"].isna()][["species","label","psa_color","size_mb"]]
if len(unc):
    print("WARNING: unclassified outgroup chromosomes:\n", unc.to_string())

# ── 4. Classify outgroup C–D categories ───────────────────────────────────────
conflicts = []
def classify_outgroup_enr(row):
    if pd.isna(row["is_macro"]): return None
    psa = row["psa_color"]
    if not row["is_macro"]:
        if psa in (CHR2_PSA, CHR13_PSA):
            conflicts.append(row)
        return "true_micro"
    if psa == CHR13_PSA: return "chr13"
    if psa == CHR2_PSA:  return "chr2"
    return "macro_rem"

df_out["enr_cat"] = df_out.apply(classify_outgroup_enr, axis=1)

if conflicts:
    print("CONFLICT — micro by karyotype but Chr2/Chr13 by psa_color:")
    for c in conflicts:
        print(f"  {c['species']} chr{c['label']} psa={c['psa_color']}")
else:
    print("Conflict check: no Chr2/Chr13 homologues found among outgroup microchromosomes. ✓")

# ── 5. Classify lacertid G–H categories ───────────────────────────────────────
def classify_lacertid_enr(psa):
    if psa == CHR18_PSA:       return "true_micro"
    if psa in MICRO_DER_COLS:  return "micro_est"
    if psa == CHR13_PSA:       return "chr13"
    if psa == CHR2_PSA:        return "chr2"
    if psa in MACRO_REM_COLS:  return "macro_rem"
    return None

df_lac["enr_cat"] = df_lac["psa_color"].map(classify_lacertid_enr)

# ── 6. Per-species reference means and enrichment ─────────────────────────────
def _ref_means(df_sp):
    return (df_sp[df_sp["enr_cat"] == "macro_rem"]
            .groupby("species")[["gc", "gene_density"]].mean())

def _add_enrichment(df_sp, ref):
    df_sp = df_sp.copy()
    df_sp["gc_enr"] = np.nan
    df_sp["gd_enr"] = np.nan
    for sp, grp in df_sp.groupby("species"):
        if sp not in ref.index: continue
        df_sp.loc[grp.index, "gc_enr"] = 100.0 * (grp["gc"] - ref.loc[sp, "gc"])
        df_sp.loc[grp.index, "gd_enr"] = grp["gene_density"] - ref.loc[sp, "gene_density"]
    return df_sp

ref_out = _ref_means(df_out)
ref_lac = _ref_means(df_lac)
df_out  = _add_enrichment(df_out, ref_out)
df_lac  = _add_enrichment(df_lac, ref_lac)

# ── 7. Sensitivity check: Chr2+Chr13 added to reference ───────────────────────
ref_lac_s = (df_lac[df_lac["enr_cat"].isin({"macro_rem", "chr2", "chr13"})]
             .groupby("species")[["gc", "gene_density"]].mean())
df_lac_s  = _add_enrichment(df_lac.copy(), ref_lac_s)
me   = df_lac[df_lac["enr_cat"] == "micro_est"]
me_s = df_lac_s[df_lac_s["enr_cat"] == "micro_est"]
print("\nSensitivity — micro-derived vs remaining macros (lacertids):")
for trait, col in [("GC enrichment (pp)", "gc_enr"), ("GD enrichment (g/Mb)", "gd_enr")]:
    d  = me[col].mean(); s = me_s[col].mean()
    print(f"  {trait}: default={d:+.3f}  with Chr2+Chr13 in ref={s:+.3f}  Δ={s-d:+.3f}")

# ── 8. Export tables ──────────────────────────────────────────────────────────
cls_all = pd.concat([
    df_out[["species","chr","label","psa_color","is_macro","enr_cat",
            "size_mb","gc","gene_density","gc_enr","gd_enr"]].assign(group="outgroup"),
    df_lac[["species","chr","label","psa_color","enr_cat",
            "size_mb","gc","gene_density","gc_enr","gd_enr"]].assign(group="lacertid", is_macro=True),
], ignore_index=True)
cls_all.to_csv(os.path.join(OUTDIR, "chr_classification.tsv"), sep="\t", index=False)

def _sp_summary(df_sp, grp_label):
    s = df_sp.dropna(subset=["enr_cat"]).groupby(["species","enr_cat"]).agg(
        n=("label","count"),
        mean_gc=("gc","mean"), mean_gd=("gene_density","mean"),
        mean_gc_enr=("gc_enr","mean"), mean_gd_enr=("gd_enr","mean"),
    ).reset_index()
    s.insert(0, "group", grp_label)
    return s

summ = pd.concat([_sp_summary(df_out, "outgroup"), _sp_summary(df_lac, "lacertid")], ignore_index=True)
summ.to_csv(os.path.join(OUTDIR, "species_summary.tsv"), sep="\t", index=False)
print("\nMean enrichment per outgroup category:")
print(summ[summ["group"]=="outgroup"].groupby("enr_cat")[["mean_gc_enr","mean_gd_enr"]].mean().round(3).to_string())
print("\nMean enrichment per lacertid category:")
print(summ[summ["group"]=="lacertid"].groupby("enr_cat")[["mean_gc_enr","mean_gd_enr"]].mean().round(3).to_string())

# ── 9. Figure layout constants ────────────────────────────────────────────────
EDGE   = "#2a2a2a"
BOX_LW = 1.1
XLABEL = r"log$_{10}$ chromosome size (Mb)"

BLOCKS = [
    ("gc",           "GC content",               "%.2f",
     "gc_enr",  "GC enrichment (pp)", "%.1f"),
    ("gene_density", "Gene density (genes/Mb)",   "%.0f",
     "gd_enr",  "Gene-density enrichment\n(genes/Mb)", "%.0f"),
]

PSA_LEGEND = [
    (MACRO_PSA,  "Macro (Chr 1,3–7,9–12)"),
    (CHR2_PSA,   "Chr 2"),  (CHR13_PSA, "Chr 13"),
    ("#FDBC6B",  "Chr 8"),  ("#BB8F6E",  "Chr 14"),
    ("#F6979A",  "Chr 15"), ("#8ED0C4",  "Chr 16"),
    (CHR17_PSA,  "Chr 17"), ("#E36EAA",  "Chr Z"),
    (CHR18_PSA,  "Chr 18 (true micro)"),
]

OUT_CATS   = ["true_micro", "chr13",    "chr2",    "macro_rem"]
OUT_LABELS = ["Retained\nmicro", "Chr 13\nhomologues", "Chr 2\nhomologues", "Remaining\nmacros"]
OUT_FILLS  = [TRUE_MICRO_FILL, CHR13_PSA, CHR2_PSA, MACRO_FILL]

LAC_CATS   = ["true_micro", "micro_est", "chr13",    "chr2",    "macro_rem"]
LAC_LABELS = ["Retained\nmicro", "Micro-\nderived", "Chr 13\nhomologues", "Chr 2\nhomologues", "Remaining\nmacros"]
LAC_FILLS  = [TRUE_MICRO_FILL, MICRO_DER_FILL, CHR13_PSA, CHR2_PSA, MACRO_FILL]

CHR_LABEL_GROUPS = {
    CHR2_PSA:  "Chr 2",
    CHR13_PSA: "Chr 13",
    "#FDBC6B":  "Chr 8",
    "#E36EAA":  "Z",
}

rng = np.random.default_rng(42)

def full_box(ax):
    for sp in ax.spines.values():
        sp.set_visible(True); sp.set_linewidth(BOX_LW); sp.set_color("#333333")

# ── 10. Enrichment panel helper ───────────────────────────────────────────────
def _enrich_panel(ax, df_sp, cats, labels, fills, enr_key, enr_label, enr_fmt,
                  sp_col_map, title=None, legend_handles=None, legend_title=None):
    ax.axhline(0, color="#888888", lw=0.9, ls="--", zorder=1)

    all_enr = df_sp[enr_key].dropna()
    if len(all_enr) == 0: return

    for xi, cat in enumerate(cats):
        sub = df_sp[df_sp["enr_cat"] == cat].dropna(subset=[enr_key])
        if len(sub) == 0: continue

        jx_ind = xi + rng.uniform(-0.22, 0.22, len(sub))
        ax.scatter(jx_ind, sub[enr_key].values,
                   color=sub["psa_color"].values,
                   s=16, alpha=0.35, zorder=2, edgecolors="none")

        sp_means = sub.groupby("species")[enr_key].mean()
        for j, (sp, val) in enumerate(sp_means.items()):
            jxi = xi + rng.uniform(-0.12, 0.12)
            ax.scatter([jxi], [val], color=sp_col_map.get(sp, "#333"),
                       s=58, edgecolors="#333", linewidths=0.45, zorder=5, alpha=0.92)

    ymin = all_enr.min(); ymax = all_enr.max()
    yr   = ymax - ymin if ymax > ymin else 1.0
    ax.set_ylim(ymin - yr*0.14, ymax + yr*0.16)
    ax.set_xticks(range(len(cats)))
    ax.set_xticklabels(labels, fontsize=8.5)
    ax.set_xlim(-0.55, len(cats) - 0.45)
    ax.set_ylabel(enr_label, fontsize=9.5)
    ax.tick_params(labelsize=8.5)
    ax.yaxis.set_major_formatter(FormatStrFormatter(enr_fmt))
    full_box(ax)
    if title:
        ax.set_title(title, fontsize=10.5, fontweight="bold", pad=6)
    if legend_handles:
        ax.legend(handles=legend_handles, title=legend_title,
                  fontsize=6.0, title_fontsize=6.5, ncol=1,
                  loc="upper left", bbox_to_anchor=(1.01, 1.0),
                  framealpha=0.9, edgecolor="#ccc",
                  borderpad=0.4, handletextpad=0.3)

# ── 11. Build figure ──────────────────────────────────────────────────────────
XMIN     = df["log_size"].min() - 0.12
XMAX     = df["log_size"].max() + 0.12
XMIN_LAC = df_lac["log_size"].min() - 0.12
XMAX_LAC = df_lac["log_size"].max() + 0.12

fig = plt.figure(figsize=(15, 15))
outer = gridspec.GridSpec(2, 3, figure=fig,
                          width_ratios=[1.45, 1.45, 1.95],
                          height_ratios=[1.0, 1.0],
                          hspace=0.26, wspace=0.33,
                          top=0.95, bottom=0.08, left=0.07, right=0.93)

ax_top    = [fig.add_subplot(outer[0, c]) for c in range(2)]
gs_out    = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[0, 2], hspace=0.52)
ax_out    = [fig.add_subplot(gs_out[i]) for i in range(2)]

ax_lac_sc = [fig.add_subplot(outer[1, c]) for c in range(2)]
gs_lac_enr = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=outer[1, 2], hspace=0.52)
ax_lac_enr = [fig.add_subplot(gs_lac_enr[i]) for i in range(2)]

# ── A–B: all-species scatter ──────────────────────────────────────────────────
for ci, (y_key, y_label, y_fmt, *_) in enumerate(BLOCKS):
    ax = ax_top[ci]
    ax.scatter(df_lac["log_size"], df_lac[y_key], s=55, color=LAC_COL,
               edgecolors="white", linewidths=0.3, alpha=1.0, zorder=2)
    for sp in ALL_OUTGROUP_SP:
        sub = df[df["species"] == sp]
        ax.scatter(sub["log_size"], sub[y_key], s=55, color=SP_COLOR[sp],
                   edgecolors="white", linewidths=0.3, alpha=0.92, zorder=3)
    ax.set_xlim(XMIN, XMAX); ax.set_ylabel(y_label, fontsize=11)
    ax.set_xlabel(XLABEL, fontsize=10.5); ax.tick_params(labelsize=9); full_box(ax)
    if y_key == "gc": ax.set_ylim(bottom=0.38)

ax_top[0].legend(
    handles=[mlines.Line2D([], [], marker="o", color="w",
                           markerfacecolor=c, markeredgecolor="none",
                           markersize=9, label=lab) for c, lab in [
        (LAC_COL,        "Lacertidae"),
        (TOXICOFERA_COL, "Toxicofera outgroups"),
        (RHINEURA_COL,   r"$\it{R.\ floridana}$"),
        (ASPIDO_COL,     r"$\it{A.\ tigris}$"),
        (GREY_MID,       r"$\it{D.\ smithi}$"),
        (GREY_DARK,      r"$\it{G.\ gallus}$"),
    ]],
    title="Lineage", title_fontsize=8.5, fontsize=7.5,
    loc="upper right", framealpha=0.95, edgecolor="#ccc")

# ── C–D: outgroup enrichment ──────────────────────────────────────────────────
out_sp_legend = [
    mlines.Line2D([], [], marker="o", color="w",
                  markerfacecolor=SP_COLOR[sp], markeredgecolor="#333",
                  linewidth=0.5, markersize=7, label=lbl)
    for sp, lbl in [
        ("Dibamus_smithi",     r"$\it{D.\ smithi}$"),
        ("Vipera_latastei",    r"$\it{V.\ latastei}$"),
        ("Iguana_iguana",      r"$\it{I.\ iguana}$"),
        ("Varanus_acanthurus", r"$\it{V.\ acanthurus}$"),
        ("Rhineura_floridana", r"$\it{R.\ floridana}$"),
        ("Aspidocelis_tigris", r"$\it{A.\ tigris}$"),
    ]
]

for ci, (_, _, _, enr_key, enr_label, enr_fmt) in enumerate(BLOCKS):
    lh  = out_sp_legend if ci == 0 else None
    ttl = "Squamate outgroups" if ci == 0 else None
    _enrich_panel(ax_out[ci], df_out, OUT_CATS, OUT_LABELS, OUT_FILLS,
                  enr_key, enr_label, enr_fmt, SP_COLOR,
                  title=ttl, legend_handles=lh, legend_title="Species")
    if ci == 1:
        ax_out[ci].text(0.01, -0.28,
            "Large points = species-level group means (remaining-macro means = 0 by construction).\n"
            "Faint points = individual chromosomes.",
            transform=ax_out[ci].transAxes,
            fontsize=6.5, color="#555", va="top")

# ── E–F: lacertid scatter with chromosome labels ──────────────────────────────
for ci, (y_key, y_label, y_fmt, *_) in enumerate(BLOCKS):
    ax = ax_lac_sc[ci]
    is_macro = df_lac["psa_color"].isin(MACRO_REM_COLS)
    for chunk, zo in [(df_lac[is_macro], 3), (df_lac[~is_macro], 4)]:
        ax.scatter(chunk["log_size"], chunk[y_key], s=90, color=chunk["psa_color"],
                   edgecolors="white", linewidths=0.3, alpha=1.0, zorder=zo)
    ax.set_xlim(XMIN_LAC, XMAX_LAC); ax.set_ylabel(y_label, fontsize=11)
    ax.set_xlabel(XLABEL, fontsize=10.5); ax.margins(y=0.11)
    ax.tick_params(labelsize=9); full_box(ax)
    ax.yaxis.set_major_formatter(FormatStrFormatter(y_fmt))

    for psa_col, lbl in CHR_LABEL_GROUPS.items():
        cluster = df_lac[df_lac["psa_color"] == psa_col].dropna(subset=[y_key])
        if len(cluster) == 0: continue
        cx = cluster["log_size"].mean()
        cy = cluster[y_key].mean()
        ax.annotate(lbl, xy=(cx, cy),
                    xytext=(6, 4), textcoords="offset points",
                    fontsize=7.5, color=psa_col, fontweight="bold",
                    arrowprops=dict(arrowstyle="-", color=psa_col, lw=0.7),
                    clip_on=True)

# ── G–H: lacertid enrichment ──────────────────────────────────────────────────
lac_sp_legend = [
    mlines.Line2D([], [], marker="o", color="w",
                  markerfacecolor=LAC_SP_COL[sp], markeredgecolor="#333",
                  linewidth=0.5, markersize=7, label=lbl)
    for sp, lbl in [
        ("Psammodromus_algirus",       r"$\it{P.\ algirus}$"),
        ("Acanthodactylus_erythrurus", r"$\it{A.\ erythrurus}$"),
        ("Eremias_argus",              r"$\it{E.\ argus}$"),
        ("Timon_lepidus",              r"$\it{T.\ lepidus}$"),
        ("Lacerta_agilis",             r"$\it{L.\ agilis}$"),
        ("Zootoca_vivipara",           r"$\it{Z.\ vivipara}$"),
        ("Scelarcis_perspicillata",    r"$\it{S.\ perspicillata}$"),
        ("Podarcis_muralis",           r"$\it{P.\ muralis}$"),
    ]
]

for ci, (_, _, _, enr_key, enr_label, enr_fmt) in enumerate(BLOCKS):
    lh  = lac_sp_legend if ci == 0 else None
    ttl = "Lacertidae" if ci == 0 else None
    _enrich_panel(ax_lac_enr[ci], df_lac, LAC_CATS, LAC_LABELS, LAC_FILLS,
                  enr_key, enr_label, enr_fmt, LAC_SP_COL,
                  title=ttl, legend_handles=lh, legend_title="Species")
    if ci == 1:
        ax_lac_enr[ci].text(0.01, -0.28,
            "Large points = species-level group means (remaining-macro means = 0 by construction).\n"
            "Faint points = individual chromosomes.",
            transform=ax_lac_enr[ci].transAxes,
            fontsize=6.5, color="#555", va="top")

# ── panel letters ─────────────────────────────────────────────────────────────
for letter, ax in zip("AB", ax_top):
    ax.text(-0.10, 1.03, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="left")
for letter, ax in zip("CD", ax_out):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="left")
for letter, ax in zip("EF", ax_lac_sc):
    ax.text(-0.10, 1.03, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="left")
for letter, ax in zip("GH", ax_lac_enr):
    ax.text(-0.14, 1.06, letter, transform=ax.transAxes,
            fontsize=14, fontweight="bold", va="bottom", ha="left")

# ── bottom synteny legend ─────────────────────────────────────────────────────
fig.legend(
    handles=[mlines.Line2D([], [], marker="o", color="w",
                           markerfacecolor=c, markeredgecolor=EDGE,
                           markersize=8, label=lab) for c, lab in PSA_LEGEND],
    title="P. algirus chromosome synteny",
    loc="lower center", ncol=5,
    fontsize=7.5, title_fontsize=8.5,
    bbox_to_anchor=(0.46, -0.025), framealpha=0.95, edgecolor="#ccc")

# ── save ──────────────────────────────────────────────────────────────────────
for ext in ("pdf", "png"):
    out = os.path.join(OUTDIR, f"plot_gc_gene_density.{ext}")
    fig.savefig(out, format=ext, bbox_inches="tight", dpi=150)
    print(f"Saved: {out}")
plt.close(fig)
