#!/usr/bin/env bash
# Run AGORA ancestral genome reconstruction.
# Requires the AGORA conda environment with agora-generic.py available.
#
# Inputs (set paths below):
#   TREE     species tree in Newick format
#   OG_DIR   directory with orthologyGroups.<NODE>.list files
#   GENES_DIR directory with genes.<species>.list files
#   WORKDIR  AGORA working/output directory
set -euo pipefail

AGORA_SCRIPT=agora-generic.py   # path to agora-generic.py
TREE=species_tree.nwk
OG_DIR=orthologyGroups
GENES_DIR=genes
WORKDIR=agora_output
THREADS=40

mkdir -p "$WORKDIR"

python "$AGORA_SCRIPT" \
    "$TREE" \
    "${OG_DIR}/orthologyGroups.%s.list" \
    "${GENES_DIR}/genes.%s.list" \
    -workingDir="$WORKDIR" \
    -nbThreads="$THREADS"
