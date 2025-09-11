#!/usr/bin/env bash

# Developed by MPA Coetzee (Forestry and Agricultural Biotechnology Institute, University of Pretoria), with the help of ChatGPT 5.

# Calculate IQ-TREE gene/site concordance factors (gCF/SCF)
# MODE options:
#   MODE=dir  : use a folder of per-gene alignments (ALN_DIR)
#   MODE=file : use one concatenated alignment (ALN_FILE) + partition file (PARTITION_FILE)
#
# Tunables via environment variables:
#   MODE=${MODE:-dir}
#   ALN_DIR=${ALN_DIR:-alignments_trimmed}
#   ALN_FILE=${ALN_FILE:-supermatrix.faa}
#   PARTITION_FILE=${PARTITION_FILE:-partitions.txt}
#   OUTDIR=${OUTDIR:-results_concordance}
#   THREADS=${THREADS:-AUTO}        # "AUTO" or an integer
#   BOOTSTRAPS=${BOOTSTRAPS:-1000}  # ultrafast bootstrap replicates
#   SCFL_SITES=${SCFL_SITES:-100}   # number of sites to sample per branch for SCF (requires >= v2.2.2)

set -euo pipefail

MODE=${MODE:-dir}
ALN_DIR=${ALN_DIR:-alignments_trimmed}
ALN_FILE=${ALN_FILE:-supermatrix.fas}
PARTITION_FILE=${PARTITION_FILE:-partitions.txt}
OUTDIR=${OUTDIR:-results_concordance}
THREADS=${THREADS:-AUTO}
BOOTSTRAPS=${BOOTSTRAPS:-1000}
SCFL_SITES=${SCFL_SITES:-100}

mkdir -p "$OUTDIR"

# Find IQ-TREE (prefer v3 if present, otherwise v2)
IQTREE_BIN=$(command -v iqtree3 || command -v iqtree2 || command -v iqtree || true)
if [[ -z "${IQTREE_BIN}" ]]; then
  echo "Error: iqtree3/iqtree2/iqtree not found in PATH." >&2
  exit 1
fi

# Normalise threads for IQ-TREE (“AUTO” or an integer)
iqtree_threads() {
  shopt -s nocasematch
  if [[ "${THREADS}" == "auto" || "${THREADS}" == "AUTO" ]]; then
    echo "AUTO"
  else
    echo "${THREADS}"
  fi
  shopt -u nocasematch
}

# Helper to test for a non-empty file
has() { [[ -s "$1" ]]; }

echo "=== IQ-TREE concordance run ==="
echo "MODE=${MODE}"
echo "OUTDIR=${OUTDIR}"
echo "Using: ${IQTREE_BIN}"
echo

case "$MODE" in
  dir)
    # Validate directory & that it contains plausible alignment files
    if [[ ! -d "$ALN_DIR" ]]; then
      echo "Error: ALN_DIR not found: $ALN_DIR" >&2
      exit 1
    fi
    shopt -s nullglob
    files=( "$ALN_DIR"/*.{fa,faa,fas,fasta,phy,phylip,aln,nex,nexus} )
    if (( ${#files[@]} == 0 )); then
      echo "Error: No alignment files detected in $ALN_DIR" >&2
      exit 1
    fi
    shopt -u nullglob

    # 1) Concatenation-based species tree (edge-linked partition model)
    if ! has "$OUTDIR/concat.treefile"; then
      echo "[1/4] Building concatenation species tree from directory: $ALN_DIR"
      "$IQTREE_BIN" -p "$ALN_DIR" \
        --prefix "$OUTDIR/concat" \
        -B "$BOOTSTRAPS" \
        -T "$(iqtree_threads)"
    else
      echo "[1/4] Skipping: $OUTDIR/concat.treefile exists."
    fi

    # 2) Per-locus trees
    if ! has "$OUTDIR/loci.treefile"; then
      echo "[2/4] Inferring locus trees from directory: $ALN_DIR"
      "$IQTREE_BIN" -S "$ALN_DIR" \
        --prefix "$OUTDIR/loci" \
        -T "$(iqtree_threads)"
    else
      echo "[2/4] Skipping: $OUTDIR/loci.treefile exists."
    fi

    # 3) Gene concordance factors
    echo "[3/4] Computing gene concordance factors (gCF)"
    "$IQTREE_BIN" -t "$OUTDIR/concat.treefile" \
      --gcf "$OUTDIR/loci.treefile" \
      --prefix "$OUTDIR/concord"

    # 4) Site concordance factors (likelihood-based; IQ-TREE >= 2.2.2)
    #echo "[4/4] Computing site concordance factors (SCF; ${SCFL_SITES} sites/branch)"
    #"$IQTREE_BIN" -te "$OUTDIR/concat.treefile" \
    #  -p "$ALN_DIR" \
    #  --scfl "$SCFL_SITES" \
    #  --prefix "$OUTDIR/concord2"
    ;;

  file)
    # Validate files
    if [[ ! -s "$ALN_FILE" ]]; then
      echo "Error: ALN_FILE not found or empty: $ALN_FILE" >&2
      exit 1
    fi
    if [[ ! -s "$PARTITION_FILE" ]]; then
      echo "Error: PARTITION_FILE not found or empty: $PARTITION_FILE" >&2
      exit 1
    fi

    # 1) Concatenation-based species tree (edge-linked partition model)
    if ! has "$OUTDIR/concat.treefile"; then
      echo "[1/4] Building concatenation species tree from file: $ALN_FILE"
      "$IQTREE_BIN" -s "$ALN_FILE" \
        -p "$PARTITION_FILE" \
        --prefix "$OUTDIR/concat" \
        -B "$BOOTSTRAPS" \
        -T "$(iqtree_threads)"
    else
      echo "[1/4] Skipping: $OUTDIR/concat.treefile exists."
    fi

    # 2) Per-locus trees
    if ! has "$OUTDIR/loci.treefile"; then
      echo "[2/4] Inferring locus trees from partition file: $PARTITION_FILE"
      "$IQTREE_BIN" -s "$ALN_FILE" \
        -S "$PARTITION_FILE" \
        --prefix "$OUTDIR/loci" \
        -T "$(iqtree_threads)"
    else
      echo "[2/4] Skipping: $OUTDIR/loci.treefile exists."
    fi

    # 3) Gene concordance factors
    echo "[3/4] Computing gene concordance factors (gCF)"
    "$IQTREE_BIN" -t "$OUTDIR/concat.treefile" \
      --gcf "$OUTDIR/loci.treefile" \
      --prefix "$OUTDIR/concord"

 #   # 4) Site concordance factors (likelihood-based; IQ-TREE >= 2.2.2)
 #   echo "[4/4] Computing site concordance factors (SCF; ${SCFL_SITES} sites/branch)"
 #   "$IQTREE_BIN" -te "$OUTDIR/concat.treefile" \
 #     -s "$ALN_FILE" \
 #     --scfl "$SCFL_SITES" \
 #     --prefix "$OUTDIR/concord2"
    ;;

  *)
    echo "Error: Unknown MODE='$MODE'. Use MODE=dir or MODE=file." >&2
    exit 1
    ;;
esac

echo
echo "Done. Key outputs in: $OUTDIR"
echo "  - concat.treefile   : concatenation-based species tree"
echo "  - loci.treefile     : set of locus trees"
echo "  - concord*          : gCF outputs/logs"
#echo "  - concord2*         : SCF outputs/logs"
