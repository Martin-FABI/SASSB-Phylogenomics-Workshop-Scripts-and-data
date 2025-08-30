#!/usr/bin/env bash

set -euo pipefail

GDIR=${GDIR:-results/gene_trees}
OUTDIR=${OUTDIR:-results/astral}

mkdir -p "$OUTDIR"

GENE_TREES="${GDIR}/gene_trees.all.tre"
BOOT_TREES="${GDIR}/gene_trees.ufboot.tre"

#Below we check for astral-pro3, you can change it to astral if you used the older version, or astral4 if you want to use the current version of astral
ASTRAL_BIN=$(command -v astral-pro3 || true)

if [[ -z "${ASTRAL_BIN}" ]]; then echo "astral not found"; exit 1; fi

#Below you can set the settings for astral
"${ASTRAL_BIN}" -i "$GENE_TREES" -o "$OUTDIR/astral.tree" 2> "$OUTDIR/astral.log"

if [[ -s "$BOOT_TREES" ]] ; then
  "${ASTRAL_BIN}" -i "$GENE_TREES" -b "$BOOT_TREES" -o "$OUTDIR/astral.mlbs.tree" 2> "$OUTDIR/astral.mlbs.log"
fi

echo "ASTRAL outputs in $OUTDIR"
