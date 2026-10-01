#!/usr/bin/env bash
set -euo pipefail

REPO="${1:-HaipeiL/prib-amp-design}"
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"

for f in checkpoint/m2_step_0010000.pt \
         checkpoint/pmic_seed_20260903_step_008000.pt \
         checkpoint/pmic_seed_20260917_step_005000.pt \
         checkpoint/pmic_seed_20261001_step_007000.pt; do
  test -s "$f" || { echo "MISSING_WEIGHT=$f"; exit 20; }
done

(cd checkpoint && sha256sum -c expected.sha256)

command -v git >/dev/null || { echo "ERROR: git missing"; exit 21; }
command -v git-lfs >/dev/null || git lfs version >/dev/null 2>&1 || { echo "ERROR: git-lfs missing"; exit 22; }
command -v gh >/dev/null || { echo "ERROR: GitHub CLI (gh) missing"; exit 23; }

gh auth status

git init
git branch -M main
git lfs install
git lfs track "checkpoint/*.pt"

git add .
git commit -m "Release PRIB AMP Design frozen research portfolio" || true

if gh repo view "$REPO" >/dev/null 2>&1; then
  git remote remove origin >/dev/null 2>&1 || true
  git remote add origin "https://github.com/${REPO}.git"
  git push -u origin main
else
  gh repo create "$REPO" \
    --public \
    --description "Adaptive antimicrobial-peptide generation with pre-experimental robustness and panel-level screening" \
    --source=. \
    --remote=origin \
    --push
fi

echo "PUBLIC_REPO=https://github.com/${REPO}"
echo "PUBLISH=PASS"
