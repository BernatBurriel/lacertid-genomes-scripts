#!/usr/bin/env bash
# Pairwise macrosynteny inference with JCVI/MCScan.
# For each species pair defined in ANCHOR_FILE, produces:
#   <query>.<subject>.blast, .anchors, .anchors.simple
#
# Inputs
#   GFF_DIR/   <species>.gff3  (longest-isoform GFF3)
#   PROT_DIR/  <species>.faa   (longest-isoform proteins)
#   SPP_FILE   one species code per line
#   ANCHOR_FILE  two columns: query  subject

GFF_DIR=gffs          # directory with <sp>.gff3 files
PROT_DIR=proteomes    # directory with <sp>.faa files
OUTDIR=macrosynteny   # output directory
SPP_FILE=species.txt
ANCHOR_FILE=anchor_pairs.txt
THREADS=32
N=10                  # parallel jobs

mkdir -p "$OUTDIR"

# 1. Convert GFF3 → BED (gene positions for JCVI)
while read -r sp; do
    [[ -f "$OUTDIR/${sp}.bed" ]] && continue
    python -m jcvi.formats.gff bed \
        --type=mRNA --key=ID \
        "$GFF_DIR/${sp}.gff3" \
        -o "$OUTDIR/${sp}.bed" &
    [[ $(jobs -r -p | wc -l) -ge $N ]] && wait -n
done < "$SPP_FILE"
wait

# 2. DIAMOND all-vs-all BLAST
while read -r query subject; do
    [[ -f "$OUTDIR/${query}.${subject}.blast" ]] && continue
    diamond makedb --in "$PROT_DIR/${subject}.faa" -d "$OUTDIR/${subject}"
    diamond blastp \
        -q "$PROT_DIR/${query}.faa" \
        -d "$OUTDIR/${subject}" \
        --outfmt 6 qseqid sseqid pident length evalue bitscore \
        --evalue 1e-5 --max-target-seqs 5 --threads "$THREADS" \
        > "$OUTDIR/${query}.${subject}.blast"
done < "$ANCHOR_FILE"

# 3. Synteny anchors (JCVI scan)
while read -r query subject; do
    [[ -f "$OUTDIR/${query}.${subject}.anchors" ]] && continue
    python -m jcvi.compara.synteny scan \
        "$OUTDIR/${query}.${subject}.blast" \
        "$OUTDIR/${query}.${subject}.anchors" \
        --qbed="$OUTDIR/${query}.bed" \
        --sbed="$OUTDIR/${subject}.bed" \
        --no_strip_names &
    [[ $(jobs -r -p | wc -l) -ge $N ]] && wait -n
done < "$ANCHOR_FILE"
wait

# 4. Screen blocks (minspan=50 as used in the paper)
while read -r query subject; do
    python -m jcvi.compara.synteny screen \
        "$OUTDIR/${query}.${subject}.anchors" \
        "$OUTDIR/${query}.${subject}.anchors.simple" \
        --qbed="$OUTDIR/${query}.bed" \
        --sbed="$OUTDIR/${subject}.bed" \
        --simple --minspan=50 &
    [[ $(jobs -r -p | wc -l) -ge $N ]] && wait -n
done < "$ANCHOR_FILE"
wait
