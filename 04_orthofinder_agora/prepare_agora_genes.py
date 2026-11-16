#!/usr/bin/env python3
"""
Convert a GFF3 annotation to the AGORA genes list format (5-column TSV):
  chromosome  start  end  orientation(+1/-1)  species|transcript_id

Usage:
    python prepare_agora_genes.py <species_name> <input.gff3[.gz]> <output.list>
"""
import sys, gzip

def opengz(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)

def parse_attrs(attr):
    return {k: v for item in attr.strip().split(";") if "=" in item
            for k, v in [item.split("=", 1)]}

if len(sys.argv) != 4:
    sys.exit("Usage: prepare_agora_genes.py <species> <input.gff3[.gz]> <output.list>")

species, gff_path, out_path = sys.argv[1], sys.argv[2], sys.argv[3]
rows = []

with opengz(gff_path) as fh:
    for line in fh:
        if not line.strip() or line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 9 or f[2] not in ("mRNA", "transcript"):
            continue
        a = parse_attrs(f[8])
        tid = a.get("ID")
        if not tid:
            continue
        orient = "1" if f[6] == "+" else "-1" if f[6] == "-" else None
        if orient is None:
            continue
        rows.append((f[0], int(f[3]), int(f[4]), orient, f"{species}|{tid}"))

rows.sort(key=lambda x: (x[0], x[1]))
with open(out_path, "w") as out:
    for r in rows:
        out.write("\t".join(map(str, r)) + "\n")

print(f"Wrote {len(rows)} transcripts to {out_path}", file=sys.stderr)
