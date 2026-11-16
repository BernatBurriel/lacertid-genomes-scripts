#!/usr/bin/env bash
# Run OrthoFinder on a directory of protein FASTA files.
# Proteins must be longest-isoform only (one sequence per gene).
# PROT_DIR should contain one <species>.faa per species.
set -euo pipefail

PROT_DIR=proteomes      # directory with <species>.faa files
OUTDIR=orthofinder      # OrthoFinder output directory
SPECIES_TREE=           # optional: path to a pre-computed species tree (.nwk)
THREADS=40
PARALLEL=16

mkdir -p "$OUTDIR"

# Copy/symlink proteomes into the OrthoFinder input directory
INPUT_DIR="${OUTDIR}/input"
mkdir -p "$INPUT_DIR"
for faa in "$PROT_DIR"/*.faa; do
    ln -sf "$(realpath "$faa")" "$INPUT_DIR/$(basename "$faa")"
done

if [[ -n "$SPECIES_TREE" ]]; then
    orthofinder \
        -f "$INPUT_DIR" \
        -S diamond \
        -M msa \
        -s "$SPECIES_TREE" \
        -t "$THREADS" \
        -a "$PARALLEL"
else
    orthofinder \
        -f "$INPUT_DIR" \
        -S diamond \
        -M msa \
        -t "$THREADS" \
        -a "$PARALLEL"
fi
