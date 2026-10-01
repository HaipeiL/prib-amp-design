#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
rm -rf generate
uv run --no-sync generate >/tmp/prib_gen1.log
sha256sum generate/library.fasta generate/top.fasta > /tmp/prib_hash1
rm -rf generate
uv run --no-sync generate >/tmp/prib_gen2.log
sha256sum generate/library.fasta generate/top.fasta > /tmp/prib_hash2
diff -u /tmp/prib_hash1 /tmp/prib_hash2
uv run --no-sync python tests/verify_frozen_submission.py
echo "REPRODUCIBILITY=PASS"
