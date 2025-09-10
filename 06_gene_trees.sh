
#!/usr/bin/env bash
set -euo pipefail

FILTER_DIR=${TRIM_DIR:-alignments_trimmed_filtered}
OUTDIR=${OUTDIR:-results/gene_trees}
THREADS=${THREADS:-AUTO}

mkdir -p "$OUTDIR"

IQTREE_BIN=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)
if [[ -z "${IQTREE_BIN}" ]]; then
  echo "iqtree3/iqtree2/iqtree not found" >&2
  exit 1
fi

# --- Cleanup: delete anything that isn't a FASTA we want (.fa .faa .fasta .fas) ---
find "$FILTER_DIR" -type f \
  \( ! -iname '*.fa' -a ! -iname '*.faa' -a ! -iname '*.fasta' -a ! -iname '*.fas' \) \
  -print -delete

# Ensure we have at least one FASTA to process
if ! find "$FILTER_DIR" -type f \( -iname '*.fa' -o -iname '*.faa' -o -iname '*.fasta' -o -iname '*.fas' \) | grep -q .; then
  echo "Error: No FASTA files (*.fa, *.faa, *.fasta, *.fas) found in $FILTER_DIR after cleanup." >&2
  exit 1
fi

shopt -s nullglob nocaseglob

# Pick the correct threads flag for your iqtree version
if "$IQTREE_BIN" -h 2>&1 | grep -qE '\-nt\b'; then
  nt_flag=(-nt "$THREADS")
else
  nt_flag=(-T "$THREADS")
fi

# Process all accepted FASTA extensions (case-insensitive)
for aln in "$FILTER_DIR"/*.fa "$FILTER_DIR"/*.faa "$FILTER_DIR"/*.fasta "$FILTER_DIR"/*.fas; do
  [[ -e "$aln" ]] || continue
  fname=$(basename "$aln")      # keep extension to avoid prefix collisions
  pre="$OUTDIR/${fname}"
  if [[ ! -s "${pre}.treefile" ]]; then
    "$IQTREE_BIN" -s "$aln" -m MFP "${nt_flag[@]}" -pre "$pre"
  fi
done

# Concatenate results only if they exist
if compgen -G "$OUTDIR/*.treefile" > /dev/null; then
  cat "$OUTDIR"/*.treefile > "$OUTDIR/gene_trees.all.tre"
fi

if compgen -G "$OUTDIR/*.ufboot" > /dev/null; then
  cat "$OUTDIR"/*.ufboot > "$OUTDIR/gene_trees.ufboot.tre"
fi

echo "TREE in $OUTDIR"
