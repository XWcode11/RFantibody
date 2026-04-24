#!/bin/bash
# End-to-end validation: RFdiffusion -> ProteinMPNN (baseline) + AntiFold (new backend)
set -e

PROJECT_ROOT="/mnt/data1/xiongw/Projects/active/AI_Protein/RFantibody"
cd "$PROJECT_ROOT"

source .venv/bin/activate
export PYTHONPATH="/tmp/antifold_repo:${PROJECT_ROOT}/src:${PROJECT_ROOT}/include/SE3Transformer:${PYTHONPATH:-}"

OUTDIR="/tmp/rfantibody_e2e"
RFD_OUT="${OUTDIR}/1_rfdiffusion"
MPNN_OUT="${OUTDIR}/2_proteinmpnn"
ANTIFOLD_OUT="${OUTDIR}/2_antifold"
mkdir -p "$RFD_OUT" "$MPNN_OUT" "$ANTIFOLD_OUT"

echo "======================================================================"
echo "[Stage 1/3] RFdiffusion: 2 antibody backbones targeting RSV site3"
echo "======================================================================"
t0=$(date +%s)
rfdiffusion \
    --target scripts/examples/example_inputs/rsv_site3.pdb \
    --framework scripts/examples/example_inputs/hu-4D5-8_Fv.pdb \
    --output "${RFD_OUT}/ab_des" \
    --num-designs 2 \
    --design-loops "L1:8-13,L2:7,L3:9-11,H1:7,H2:6,H3:5-13" \
    --hotspots "T305,T456" \
    --diffuser-t 50 \
    --deterministic 2>&1 | tail -30
echo "[Stage 1 done in $(($(date +%s)-t0))s]"
echo ""
echo "RFdiffusion outputs:"
ls -l "${RFD_OUT}"/*.pdb 2>/dev/null | head -10

echo ""
echo "======================================================================"
echo "[Stage 2a/3] ProteinMPNN (baseline backend)"
echo "======================================================================"
t0=$(date +%s)
proteinmpnn \
    --input-dir "${RFD_OUT}" \
    --output-dir "${MPNN_OUT}" \
    --seqs-per-struct 2 \
    --temperature 0.2 \
    --backend proteinmpnn 2>&1 | tail -20
echo "[Stage 2a done in $(($(date +%s)-t0))s]"
echo ""
echo "ProteinMPNN outputs:"
ls -l "${MPNN_OUT}"/*.pdb 2>/dev/null | head -10

echo ""
echo "======================================================================"
echo "[Stage 2b/3] AntiFold (new backend)"
echo "======================================================================"
t0=$(date +%s)
proteinmpnn \
    --input-dir "${RFD_OUT}" \
    --output-dir "${ANTIFOLD_OUT}" \
    --seqs-per-struct 2 \
    --temperature 0.2 \
    --backend antifold 2>&1 | tail -20
echo "[Stage 2b done in $(($(date +%s)-t0))s]"
echo ""
echo "AntiFold outputs:"
ls -l "${ANTIFOLD_OUT}"/*.pdb 2>/dev/null | head -10

echo ""
echo "======================================================================"
echo "[Summary] H-chain sequence per design (truncated)"
echo "======================================================================"
for pdb in "${RFD_OUT}"/*.pdb; do
    [ -f "$pdb" ] || continue
    tag=$(basename "$pdb" .pdb)
    echo ""
    echo "### ${tag}"
    for backend_dir in "${MPNN_OUT}" "${ANTIFOLD_OUT}"; do
        bname=$(basename "$backend_dir")
        for out_pdb in "${backend_dir}/${tag}"_dldesign_*.pdb; do
            [ -f "$out_pdb" ] || continue
            seq=$(python3 -c "
aa3_to_1={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
seen=set(); seq=''
with open('$out_pdb') as f:
    for line in f:
        if line.startswith('ATOM') and line[12:16].strip()=='CA' and line[21]=='H':
            k=(line[21],line[22:27])
            if k not in seen: seen.add(k); seq+=aa3_to_1.get(line[17:20].strip(),'X')
print(seq)
")
            printf "  [%-15s] %s: %s\n" "$bname" "$(basename "$out_pdb" .pdb)" "$seq"
        done
    done
done

echo ""
echo "======================================================================"
echo "E2E pipeline validation complete."
echo "  RFdiffusion backbones: ${RFD_OUT}"
echo "  ProteinMPNN designs:   ${MPNN_OUT}"
echo "  AntiFold designs:      ${ANTIFOLD_OUT}"
echo "======================================================================"
