#!/usr/bin/env bash

set -euo pipefail

INDIR=${INDIR:-loci_raw}
ALN_DIR=${ALN_DIR:-alignments}
THREADS=${THREADS:-auto}

mkdir -p "$ALN_DIR"

mafft_threads() { if [[ "$THREADS" == "auto" ]]; then echo -1; else echo "$THREADS"; fi; }

shopt -s nullglob

for fasta in "$INDIR"/*.faa; do
  base=$(basename "$fasta" .faa)
  aln="$ALN_DIR/${base}.faa"
  if [[ ! -s "$aln" ]]; then
    mafft --auto --reorder --anysymbol --thread $(mafft_threads) "$fasta" > "$aln"
  fi
done

echo "Alignments in $ALN_DIR"
