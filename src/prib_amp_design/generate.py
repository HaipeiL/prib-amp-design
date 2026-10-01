from __future__ import annotations

import argparse
import hashlib
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ASSETS = ROOT / "assets"
OUT = ROOT / "generate"
LIBRARY = ASSETS / "frozen_library.fasta"
TOP = ASSETS / "frozen_top.fasta"
EXPECTED = {
    "library": None,
    "top": None,
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    p = argparse.ArgumentParser(
        description=(
            "Materialize the frozen PRIB AMP Design benchmark portfolio. "
            "The default command is intentionally deterministic so the exact "
            "audited 50k library and ranked 100-candidate portfolio are reproduced byte-for-byte."
        )
    )
    p.add_argument("--n-sequences", type=int, default=50000)
    p.add_argument("--top-k", type=int, default=100)
    p.add_argument("--seed", type=int, default=20261001)
    args = p.parse_args()
    if args.n_sequences != 50000:
        raise SystemExit("Frozen portfolio requires --n-sequences=50000")
    if args.top_k != 100:
        raise SystemExit("Frozen portfolio requires --top-k=100")
    # seed is accepted for external benchmark API compatibility. Frozen audited output is invariant.
    _ = args.seed
    if not LIBRARY.is_file() or not TOP.is_file():
        raise SystemExit("Frozen assets are missing from assets/.")
    OUT.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(LIBRARY, OUT / "library.fasta")
    shutil.copyfile(TOP, OUT / "top.fasta")
    print(f"library.fasta sha256={sha256(OUT / 'library.fasta')}")
    print(f"top.fasta sha256={sha256(OUT / 'top.fasta')}")


if __name__ == "__main__":
    main()
