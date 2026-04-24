#!/bin/bash

# Run RF2 structure prediction on ProteinMPNN and AntiFold designs for comparison.
#
# Consumes outputs produced by e2e_antifold_validation.sh (under the same
# OUTDIR) and writes per-design confidence scores for side-by-side review.

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EXAMPLES_DIR="$SCRIPT_DIR"

# Use fewer recycles for faster throughput during comparison
NUM_RECYCLES=${NUM_RECYCLES:-1}
SEED=${SEED:-42}

OUTDIR="${OUTDIR:-$EXAMPLES_DIR/example_outputs/e2e_antifold_validation}"
MPNN_IN="${OUTDIR}/2_proteinmpnn"
ANTIFOLD_IN="${OUTDIR}/2_antifold"

MPNN_OUT="${OUTDIR}/3_rf2_mpnn"
ANTIFOLD_OUT="${OUTDIR}/3_rf2_antifold"
mkdir -p "$MPNN_OUT" "$ANTIFOLD_OUT"

echo "======================================================================"
echo "[RF2] Predicting ProteinMPNN designs  -> ${MPNN_OUT}"
echo "        ${NUM_RECYCLES} recycles, seed=${SEED}"
echo "======================================================================"
uv run rf2 -i "$MPNN_IN" -o "$MPNN_OUT" -r "$NUM_RECYCLES" --seed "$SEED" --no-cautious
echo "[MPNN RF2 done]"

echo ""
echo "======================================================================"
echo "[RF2] Predicting AntiFold designs       -> ${ANTIFOLD_OUT}"
echo "        ${NUM_RECYCLES} recycles, seed=${SEED}"
echo "======================================================================"
uv run rf2 -i "$ANTIFOLD_IN" -o "$ANTIFOLD_OUT" -r "$NUM_RECYCLES" --seed "$SEED" --no-cautious
echo "[AntiFold RF2 done]"

echo ""
echo "======================================================================"
echo "[Summary] RF2 confidence scores per design"
echo "======================================================================"
for dir in "$MPNN_OUT" "$ANTIFOLD_OUT"; do
    bname=$(basename "$dir")
    echo ""
    echo "--- $bname ---"
    for pdb in "$dir"/*best.pdb; do
        [ -f "$pdb" ] || continue
        tag=$(basename "$pdb" _best.pdb)
        # Extract SCORE lines
        plddt=$(grep "SCORE pred_lddt" "$pdb" | awk '{print $NF}')
        ipae=$(grep "SCORE interaction_pae" "$pdb" | awk '{print $NF}')
        printf "  %-40s  pLDDT=%s  ipAE=%s\n" "$tag" "$plddt" "$ipae"
    done
done

echo ""
echo "======================================================================"
echo "RF2 scoring complete."
echo "  ProteinMPNN predictions: ${MPNN_OUT}"
echo "  AntiFold predictions:    ${ANTIFOLD_OUT}"
echo "======================================================================"
