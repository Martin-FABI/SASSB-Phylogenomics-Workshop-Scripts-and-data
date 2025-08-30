#!/usr/bin/env bash
set -euo pipefail

BUSCO_DIR=${BUSCO_DIR:-busco_data}
OUTDIR=${OUTDIR:-loci_raw}

# Make sure the output directory exists
mkdir -p "$OUTDIR"

# Expand globs safely (no literal patterns if nothing matches)
shopt -s nullglob

# Remember where we started so we can write outputs outside BUSCO_DIR
ROOT_DIR=$(pwd)

# Change into BUSCO_DIR so the multi-glob works even if BUSCO_DIR has spaces
pushd "$BUSCO_DIR" >/dev/null

# Collect all FASTA-like files in one pass
files=( *.fasta *.fa *.fas *.faa )

if (( ${#files[@]} == 0 )); then
  echo "No FASTA files (.fasta .fa .fas .faa) found in $BUSCO_DIR"
  popd >/dev/null
  exit 0
fi

for file in "${files[@]}"; do
  locus="${file%.*}"                         # strip last extension
  out="$ROOT_DIR/$OUTDIR/${locus}.faa"

  # Keep only text up to the FIRST '|' on header lines; leave sequences untouched
  awk '/^>/{sub(/\|.*/,"")} {print}' "$file" > "$out"
  # sed '/^>/{s/\|.*//}' "$file" > "$out"
done

popd >/dev/null

echo "Headers truncated at first '|' — outputs written to $OUTDIR"

