#!/usr/bin/env bash
#Developed by Martin Coetzee (University of Pretoria, Forestry and Acgricultural Biotechnology Insitute), with assistance from ChatGPT.

set -euo pipefail

INDIR=${INDIR:-alignments}
TRIM_DIR=${TRIM_DIR:-alignments_trimmed}
THREADS=${THREADS:-auto}

mkdir -p "$TRIM_DIR"

shopt -s nullglob

for fasta in "$INDIR"/*.faa; do
  base=$(basename "$fasta" .faa)
  aln="$INDIR/${base}.faa"
  trim="$TRIM_DIR/${base}.trim.faa"
  #  if [[ ! -s "$aln" ]]; then
  #  mafft --auto --reorder --anysymbol --thread $(mafft_threads) "$fasta" > "$aln"
  #fi
	if [[ ! -s "$trim" ]]; then
		trimal -in "$aln" -out "$trim" -automated1
	fi
done

echo "Trimmed in $TRIM_DIR"
