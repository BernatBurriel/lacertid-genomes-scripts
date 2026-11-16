#!/usr/bin/env python3
"""
Export outgroup (Aspidocelis tigris, Rhineura floridana) Hi-C contact data.

Outputs (in OUTDIR):
  hic_outgroup_trans_fractions.tsv
  hic_outgroup_log2oe_pairs.tsv
  hic_outgroup_chrom_info.tsv

Requires: cooler, numpy
"""
import numpy as np
import cooler
import os

OUTDIR   = "hic_contacts"
COOL_DIR = "cool_files"
os.makedirs(OUTDIR, exist_ok=True)

CONFIGS = {
    "Aspidocelis_tigris": {
        "cool":   f"{COOL_DIR}/Aspidocelis_tigris.50000.cool",
        "label":  "A. tigris",
        "prefix": "CM",
        "macro":  [f"CM{107880 + i}.1" for i in range(11)],    # chr 1–11
        "micro":  [f"CM{107891 + i}.1" for i in range(12)],    # chr 12–23
        "num":    {f"CM{107880 + i}.1": str(i + 1) for i in range(23)},
        "exclude": [],
    },
    "Rhineura_floridana": {
        "cool":   f"{COOL_DIR}/Rhineura_floridana.50000.cool",
        "label":  "R. floridana",
        "prefix": "Chromosome",
        "macro":  [f"Chromosome_{i}" for i in range(1, 12)],   # chr 1–11
        "micro":  [f"Chromosome_{i}" for i in range(12, 23)],  # chr 12–22
        "num":    {f"Chromosome_{i}": str(i) for i in range(1, 23)},
        "exclude": [],
    },
}


def classify(chrom, cfg):
    if chrom in cfg["micro"]:  return "micro"
    if chrom in cfg["macro"]:  return "macro"
    return None


info_rows, trans_rows, oe_rows = [], [], []

for sp, cfg in CONFIGS.items():
    if not os.path.exists(cfg["cool"]):
        print(f"Skipping {sp}: {cfg['cool']} not found", flush=True)
        continue

    print(f"Processing {sp} ...", flush=True)
    clr    = cooler.Cooler(cfg["cool"])
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

    csizes = dict(zip(clr.chromnames, clr.chromsizes))
    chroms = [c for c in clr.chromnames
              if c.startswith(cfg["prefix"])
              and c not in cfg["exclude"]
              and classify(c, cfg) is not None]

    chrom_total = {}
    for ch in chroms:
        m = (pixels["c1"] == ch) | (pixels["c2"] == ch)
        chrom_total[ch] = pixels.loc[m, "balanced"].sum()
    grand_total = sum(chrom_total.values()) / 2

    for ch in chroms:
        cls = classify(ch, cfg)
        num = cfg["num"].get(ch, ch)
        m_tot = (pixels["c1"] == ch) | (pixels["c2"] == ch)
        m_cis = (pixels["c1"] == ch) & (pixels["c2"] == ch)
        total = pixels.loc[m_tot, "balanced"].sum()
        cis   = pixels.loc[m_cis, "balanced"].sum()
        frac  = (total - cis) / total if total > 0 else float("nan")
        info_rows.append((sp, cfg["label"], ch, num, cls, csizes[ch] / 1e6))
        trans_rows.append((sp, cfg["label"], ch, num, cls, csizes[ch] / 1e6, frac))

    idx = {c: i for i, c in enumerate(chroms)}
    mat_obs = np.zeros((len(chroms), len(chroms)))
    for p in pixels.itertuples(index=False):
        c1, c2 = p.c1, p.c2
        if c1 == c2 or c1 not in idx or c2 not in idx:
            continue
        i, j = idx[c1], idx[c2]
        mat_obs[i, j] += p.balanced
        mat_obs[j, i] += p.balanced

    for i, c1 in enumerate(chroms):
        for j, c2 in enumerate(chroms):
            if j <= i:
                continue
            obs    = mat_obs[i, j]
            exp    = (chrom_total[c1] * chrom_total[c2]) / grand_total
            log2oe = np.log2(obs / exp) if obs > 0 else float("nan")
            n1, n2   = cfg["num"].get(c1, c1), cfg["num"].get(c2, c2)
            cl1, cl2 = classify(c1, cfg), classify(c2, cfg)
            oe_rows.append((sp, cfg["label"], c1, n1, cl1, c2, n2, cl2, log2oe))

with open(f"{OUTDIR}/hic_outgroup_chrom_info.tsv", "w") as f:
    f.write("species\tspecies_label\tchrom\tchrom_num\tclass\tsize_mb\n")
    for r in info_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

with open(f"{OUTDIR}/hic_outgroup_trans_fractions.tsv", "w") as f:
    f.write("species\tspecies_label\tchrom\tchrom_num\tclass\tsize_mb\ttrans_frac\n")
    for r in trans_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

with open(f"{OUTDIR}/hic_outgroup_log2oe_pairs.tsv", "w") as f:
    f.write("species\tspecies_label\tchrom1\tchrom1_num\tclass1\tchrom2\tchrom2_num\tclass2\tlog2_oe\n")
    for r in oe_rows:
        f.write("\t".join(str(x) for x in r) + "\n")

print(f"Done. Output in {OUTDIR}/")
