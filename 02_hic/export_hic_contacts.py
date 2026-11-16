#!/usr/bin/env python3
"""
Export lacertid Hi-C interchromosomal contact data from .cool files.

Outputs (in OUTDIR):
  hic_trans_fractions.tsv   per-chromosome trans-contact fraction
  hic_log2oe_pairs.tsv      log2(observed/expected) for every chromosome pair
  hic_chrom_info.tsv        chromosome metadata

Requires: cooler, numpy
"""
import numpy as np
import cooler

OUTDIR = "hic_contacts"   # output directory
COOL_DIR = "cool_files"   # directory containing .cool files

# Chromosome classification per species
CLASSIFICATIONS = {
    "Psammodromus_algirus": {
        "micro_derived": ["SUPER_8", "SUPER_14", "SUPER_15", "SUPER_16"],
        "true_micro":    ["SUPER_18"],
        "exclude":       ["SUPER_W", "SUPER_Z"],
    },
    "Acanthodactylus_erythrurus": {
        "micro_derived": ["SUPER_8", "SUPER_13", "SUPER_15", "SUPER_17"],
        "true_micro":    ["SUPER_18"],
        "exclude":       ["SUPER_Z"],
    },
    "Timon_lepidus": {
        "micro_derived": ["SUPER_7", "SUPER_12", "SUPER_14", "SUPER_15"],
        "true_micro":    ["SUPER_17"],
        "exclude":       ["SUPER_Z"],
    },
    "Scelarcis_perspicillata": {
        "micro_derived": ["SUPER_8", "SUPER_13", "SUPER_15", "SUPER_16"],
        "true_micro":    ["SUPER_18"],
        "exclude":       ["SUPER_Z"],
    },
}

COOL_PATHS = {sp: f"{COOL_DIR}/{sp}.50000.cool" for sp in CLASSIFICATIONS}

SP_LABELS = {
    "Psammodromus_algirus":    "P. algirus",
    "Acanthodactylus_erythrurus": "A. erythrurus",
    "Timon_lepidus":           "T. lepidus",
    "Scelarcis_perspicillata": "S. perspicillata",
}

import os
os.makedirs(OUTDIR, exist_ok=True)


def classify(chrom, sp_cls):
    if chrom in sp_cls.get("exclude", []):  return "excluded"
    if chrom in sp_cls["micro_derived"]:    return "micro_derived"
    if chrom in sp_cls["true_micro"]:       return "true_micro"
    return "macro_derived"


trans_rows, oe_rows, info_rows = [], [], []

for sp, cool_path in COOL_PATHS.items():
    print(f"Processing {sp} ...", flush=True)
    clr    = cooler.Cooler(cool_path)
    sp_cls = CLASSIFICATIONS[sp]
    bins   = clr.bins()[:]
    pixels = clr.pixels()[:]

    w = bins["weight"].values
    pixels["balanced"] = (pixels["count"].values
                          * w[pixels["bin1_id"].values]
                          * w[pixels["bin2_id"].values])
    pixels = pixels.dropna(subset=["balanced"])

    cob = bins["chrom"].values
    pixels["c1"] = cob[pixels["bin1_id"].values]
    pixels["c2"] = cob[pixels["bin2_id"].values]

    csizes       = dict(zip(clr.chromnames, clr.chromsizes))
    super_chroms = [c for c in clr.chromnames
                    if c.startswith("SUPER") and classify(c, sp_cls) != "excluded"]

    chrom_total = {}
    for ch in super_chroms:
        m = (pixels["c1"] == ch) | (pixels["c2"] == ch)
        chrom_total[ch] = pixels.loc[m, "balanced"].sum()
    grand_total = sum(chrom_total.values()) / 2

    for ch in super_chroms:
        m_tot = (pixels["c1"] == ch) | (pixels["c2"] == ch)
        m_cis = (pixels["c1"] == ch) & (pixels["c2"] == ch)
        total = pixels.loc[m_tot, "balanced"].sum()
        cis   = pixels.loc[m_cis, "balanced"].sum()
        frac  = (total - cis) / total if total > 0 else float("nan")
        cls   = classify(ch, sp_cls)
        trans_rows.append((sp, SP_LABELS[sp], ch, cls, csizes[ch] / 1e6, frac))
        info_rows.append((sp, SP_LABELS[sp], ch, cls, csizes[ch] / 1e6))

    idx = {c: i for i, c in enumerate(super_chroms)}
    n   = len(super_chroms)
    mat_obs = np.zeros((n, n))
    for p in pixels.itertuples(index=False):
        c1, c2 = p.c1, p.c2
        if c1 == c2 or c1 not in idx or c2 not in idx:
            continue
        i, j = idx[c1], idx[c2]
        mat_obs[i, j] += p.balanced
        mat_obs[j, i] += p.balanced

    for i, c1 in enumerate(super_chroms):
        for j, c2 in enumerate(super_chroms):
            if j <= i:
                continue
            obs    = mat_obs[i, j]
            exp    = (chrom_total[c1] * chrom_total[c2]) / grand_total
            log2oe = np.log2(obs / exp) if obs > 0 else float("nan")
            oe_rows.append((sp, SP_LABELS[sp], c1, classify(c1, sp_cls),
                            c2, classify(c2, sp_cls), log2oe))

with open(f"{OUTDIR}/hic_trans_fractions.tsv", "w") as f:
    f.write("species\tspecies_label\tchrom\tclass\tsize_mb\ttrans_frac\n")
    for r in trans_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

with open(f"{OUTDIR}/hic_log2oe_pairs.tsv", "w") as f:
    f.write("species\tspecies_label\tchrom1\tclass1\tchrom2\tclass2\tlog2_oe\n")
    for r in oe_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

with open(f"{OUTDIR}/hic_chrom_info.tsv", "w") as f:
    f.write("species\tspecies_label\tchrom\tclass\tsize_mb\n")
    for r in info_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

print(f"Done. Output in {OUTDIR}/")
