#!/usr/bin/env bash
set -euo pipefail

INDIR=${INDIR:-alignments_trimmed}
FILTER_DIR=${FILTER_DIR:-alignments_trimmed_filtered}
MINLEN=${MINLEN:-200}

mkdir -p "$FILTER_DIR"
shopt -s nullglob

command -v bioawk >/dev/null || { echo "Error: bioawk not found in PATH" >&2; exit 1; }

for fasta in "$INDIR"/*.faa; do
  # Get alignment length only if all sequences have the same length
  if len=$(bioawk -c fastx '
      NR==1 {L=length($seq)}
      { if (length($seq)!=L) exit 1 }
      END { if (NR>0) print L }' "$fasta"); then
    # Echo the alignment length once (before filtering)
    echo "$(basename "$fasta"): alignment length = $len"

    # Filter by minimum alignment length
    if (( len > MINLEN )); then
      cp -p -- "$fasta" "$FILTER_DIR"/
    fi
  fi
done

echo "Filtered in $FILTER_DIR"