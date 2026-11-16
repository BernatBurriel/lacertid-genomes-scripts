#!/usr/bin/env python3
"""
Count synteny blocks per (species, focal_chr, outgroup_chr) from JCVI .anchors.simple files.

Usage:
    python parse_synteny_blocks.py <minspan_suffix>
    e.g.  python parse_synteny_blocks.py 50   # reads .anchors.simple_ms50

Output: TSV  species | focal_chr | outgroup_chr | outgroup_class | n_blocks
"""
import os, re, sys
from collections import defaultdict

# ── Configuration ────────────────────────────────────────────────────────────
OUTDIR   = "macrosynteny"       # directory with .anchors.simple_ms* files
FOCAL_SP = "Psammodromus_algirus"
FOCAL_CHRS = {f"SUPER_{i}" for i in range(1, 19)} | {"SUPER_Z"}

# Macrochromosomes per outgroup species (all others are micro)
MACRO_CHRS = {
    "Gallus_gallus":         {f"Chromosome_{i}" for i in range(1, 10)} | {"Chromosome_Z"},
    "Dibamus_smithi":        {f"Chromosome_{i}" for i in range(1, 12)} | {"Chromosome_X", "Chromosome_Y"},
    "Varanus_acanthurus":    {f"Chromosome_{i}" for i in range(1, 9)},
    "Iguana_iguana":         {f"Chromosome_{i}" for i in range(1, 7)},
    "Vipera_latastei":       {f"SUPER_{i}" for i in range(1, 8)} | {"SUPER_Z"},
    "Aspidocelis_tigris":    {f"CM{107880 + i}.1" for i in range(11)},
    "Rhineura_floridana":    {f"Chromosome_{i}" for i in range(1, 12)},
}

SPECIES = [
    "Gallus_gallus", "Dibamus_smithi", "Varanus_acanthurus",
    "Iguana_iguana", "Vipera_latastei", "Aspidocelis_tigris", "Rhineura_floridana",
]

SCAFFOLD_PREFIXES = ("Scaffold_", "scaffold_", "ragtag_")

# ── Helpers ──────────────────────────────────────────────────────────────────
def is_scaffold(c):
    return any(c.startswith(p) for p in SCAFFOLD_PREFIXES)

def classify(sp, chrom):
    if is_scaffold(chrom):
        return None
    return "macro" if chrom in MACRO_CHRS[sp] else "micro"

def gene_to_chrom(gene):
    return gene.rsplit("-g", 1)[0]

def simplify(chrom):
    m = re.match(r"^CM(\d+)\.1$", chrom)
    if m:
        return str(int(m.group(1)) - 107879)
    return re.sub(r"^(Chromosome_|SUPER_)", "", chrom)

# ── Main ─────────────────────────────────────────────────────────────────────
suffix = sys.argv[1] if len(sys.argv) > 1 else "50"

rows = []
for sp in SPECIES:
    simple = os.path.join(OUTDIR, f"{sp}.{FOCAL_SP}.anchors.simple_ms{suffix}")
    if not os.path.exists(simple):
        print(f"[WARN] missing: {simple}", file=sys.stderr)
        continue

    counts = defaultdict(int)
    with open(simple) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            q_gene, s_gene = parts[0], parts[1]
            q_chr = gene_to_chrom(q_gene)   # outgroup
            s_chr = gene_to_chrom(s_gene)   # focal (Psam)
            if s_chr not in FOCAL_CHRS:
                continue
            cls = classify(sp, q_chr)
            if cls is None:
                continue
            counts[(s_chr, q_chr, cls)] += 1

    for (focal_chr, out_chr, cls), n in sorted(counts.items()):
        rows.append((sp, simplify(focal_chr), simplify(out_chr), cls, n))

print("species\tfocal_chr\toutgroup_chr\toutgroup_class\tn_blocks")
for r in rows:
    print("\t".join(str(x) for x in r))
