#!/usr/bin/env bash
# Sensitivity analysis: re-screen synteny blocks at minspan=30 and minspan=50.
# Produces .anchors.simple_ms30 and .anchors.simple_ms50 for each species pair.
# jcvi always writes the block summary to .anchors.simple regardless of the
# output argument; this script copies that file after each run.
set -euo pipefail

OUTDIR=macrosynteny
ANCHOR_FILE=anchor_pairs.txt

for MS in 30 50; do
    echo "[$(date +%T)] minspan=${MS}"
    while read -r query subject; do
        ANCHORS="$OUTDIR/${query}.${subject}.anchors"
        SIMPLE="$OUTDIR/${query}.${subject}.anchors.simple"
        THROWAWAY="${ANCHORS}.tmp"

        python -m jcvi.compara.synteny screen \
            "$ANCHORS" "$THROWAWAY" \
            --qbed="$OUTDIR/${query}.bed" \
            --sbed="$OUTDIR/${subject}.bed" \
            --simple --minspan="$MS" 2>/dev/null

        cp "$SIMPLE" "${SIMPLE}_ms${MS}"
        rm -f "$THROWAWAY"
        echo "  ${query} vs ${subject}: $(wc -l < ${SIMPLE}_ms${MS}) blocks"
    done < "$ANCHOR_FILE"
done
echo "Done."
