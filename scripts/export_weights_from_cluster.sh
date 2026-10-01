#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRIB_ROOT="${PRIB_ROOT:-/pfss/mlde/workspaces/mlde_wsp_PRIB_KI/PRIB}"
OUT="$ROOT/checkpoint"
mkdir -p "$OUT"
cp "$PRIB_ROOT/projects/amp_adaptive_cardinality/runs/paper_finetune_stable_v3/M2_full/checkpoints/step_0010000.pt" "$OUT/m2_step_0010000.pt"
cp "$PRIB_ROOT/projects/amp_adaptive_cardinality/runs/independent_pmic_evaluator_v1/seed_20260903/checkpoints/step_008000.pt" "$OUT/pmic_seed_20260903_step_008000.pt"
cp "$PRIB_ROOT/projects/amp_adaptive_cardinality/runs/independent_pmic_evaluator_v1/seed_20260917/checkpoints/step_005000.pt" "$OUT/pmic_seed_20260917_step_005000.pt"
cp "$PRIB_ROOT/projects/amp_adaptive_cardinality/runs/independent_pmic_evaluator_v1/seed_20261001/checkpoints/step_007000.pt" "$OUT/pmic_seed_20261001_step_007000.pt"
cat > "$OUT/expected.sha256" <<'HASHES'
fdc17a1052dcd00678bbe317dbd7e0d29ee36ba36e5cc6221707f07ffb8e2141  m2_step_0010000.pt
a690836e8f9fe68e318a012fa35bb1c80615d608ea1173466ac6e49479eb181a  pmic_seed_20260903_step_008000.pt
5eaf2024023c53fef9d5ba58ce3d375fc448a3eaad62cf3365ba6dd8f2ac1ea3  pmic_seed_20260917_step_005000.pt
2881231b3a631ba874f0a9d6a287b864f9c99a2d882130a8238878cdf58e223b  pmic_seed_20261001_step_007000.pt
HASHES
(cd "$OUT" && sha256sum -c expected.sha256)
echo "WEIGHTS_EXPORT=PASS"
