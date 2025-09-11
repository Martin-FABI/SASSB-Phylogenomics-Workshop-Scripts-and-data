#!/usr/bin/env bash

set -euo pipefail

BUSCO_DIR=${BUSCO_DIR:-busco_runs}

OUTDIR=${OUTDIR:-loci_raw}
mkdir -p "$OUTDIR"
rm -f "$OUTDIR"/*.faa || true
shopt -s nullglob

for taxon_dir in "$BUSCO_DIR"/*; do
  [[ -d "$taxon_dir" ]] || continue
  taxon=$(basename "$taxon_dir")
  for faa in "$taxon_dir"/run_*/busco_sequences/single_copy_busco/*.faa; do
    [[ -f "$faa" ]] || continue
    locus=$(basename "$faa" .faa)
    out="$OUTDIR/${locus}.faa"
    awk -v sp="$taxon" -v id="$locus" 'BEGIN{ORS="\n"} /^>/{print ">" sp "|" id; next} {print}' "$faa" >> "$out"
  done
done

echo "Assembled per-locus FASTAs in $OUTDIR"
