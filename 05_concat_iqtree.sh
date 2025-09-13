#!/usr/bin/env bash
#Developed by Martin Coetzee (University of Pretoria, Forestry and Acgricultural Biotechnology Insitute), with assistance from ChatGPT.

set -euo pipefail

# Input directory with alignments (FASTA files only after cleanup)
FILTER_DIR=${TRIM_DIR:-alignments_trimmed_filtered}
OUTROOT=${OUTROOT:-results/concat/concat}
THREADS=${THREADS:-AUTO}

mkdir -p "$(dirname "$OUTROOT")"

# Find iqtree2 or iqtree
IQTREE_BIN=$(command -v iqtree2 || command -v iqtree || true)
if [[ -z "${IQTREE_BIN}" ]]; then
  echo "iqtree2/iqtree not found" >&2
  exit 1
fi

# --- Cleanup: remove non-FASTA files that can break -S directory scans (e.g. .DS_Store) ---
# Allowed extensions: .fa .faa .fasta .fas  (case-insensitive)
# Prints what it removes, then deletes it.
find "$FILTER_DIR" -type f \
  \( ! -iname '*.fa' -a ! -iname '*.faa' -a ! -iname '*.fasta' -a ! -iname '*.fas' \) \
  -print -delete

# Sanity check: ensure we still have FASTA files to process
if ! find "$FILTER_DIR" -type f \( -iname '*.fa' -o -iname '*.faa' -o -iname '*.fasta' -o -iname '*.fas' \) | grep -q .; then
  echo "Error: No FASTA files found in $FILTER_DIR after cleanup." >&2
  exit 1
fi

# Run IQ-TREE on the directory of alignments
"${IQTREE_BIN}" -p "$FILTER_DIR" -m MFP -bb 1000 -alrt 1000 -T "${THREADS}" -pre "$OUTROOT"
#For a quick check, use below
#"${IQTREE_BIN}" -S "$FILTER_DIR" -m MFP -T "${THREADS}" -pre "$OUTROOT"

echo "TREE in $OUTROOT"
