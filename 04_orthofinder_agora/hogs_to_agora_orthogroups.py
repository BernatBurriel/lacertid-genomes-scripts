#!/usr/bin/env python3
"""
Convert OrthoFinder HOG tables (N*.tsv from Phylogenetic_Hierarchical_Orthogroups/)
to AGORA orthologyGroups format.

Usage:
    python hogs_to_agora_orthogroups.py <species_list.txt> <hog_dir> <outdir>

species_list.txt: one species name per line, matching column headers in HOG tables.
hog_dir:          directory containing N0.tsv, N1.tsv, … from OrthoFinder.
outdir:           output directory; writes orthologyGroups.<NODE>.list per node.
"""
import os, sys, csv, glob

if len(sys.argv) != 4:
    sys.exit("Usage: hogs_to_agora_orthogroups.py <species_list.txt> <hog_dir> <outdir>")

species_file, hog_dir, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
os.makedirs(outdir, exist_ok=True)

META_COLS = {"HOG", "OG", "Orthogroup", "Gene Tree Parent Clade",
             "Taxon", "Node", "Clade", "Description"}

with open(species_file) as fh:
    valid_species = {x.strip() for x in fh if x.strip()}

for tsv in sorted(glob.glob(os.path.join(hog_dir, "N*.tsv"))):
    anc = os.path.basename(tsv).replace(".tsv", "")
    out_path = os.path.join(outdir, f"orthologyGroups.{anc}.list")

    with open(tsv, newline="") as fh, open(out_path, "w") as out:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            continue
        headers = [h.strip() for h in reader.fieldnames]
        sp_cols = [h for h in headers if h not in META_COLS and h in valid_species]

        for row in reader:
            hog_id = row.get("HOG", "")
            genes  = []
            for sp in sp_cols:
                cell = row.get(sp, "").strip()
                if not cell:
                    continue
                for g in cell.split(", "):
                    g = g.strip()
                    if g:
                        genes.append(f"{sp}|{g}")
            if len(genes) >= 2:
                out.write(hog_id + "\t" + "\t".join(genes) + "\n")

    print(f"  {anc}: written {out_path}", file=sys.stderr)
