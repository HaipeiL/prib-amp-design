# External evaluation interface

PRIB AMP Design exposes a deterministic external-evaluation interface:

```bash
uv sync
uv run generate
```

The command writes an exactly reproducible 50,000-member library and ranked 100-candidate portfolio to `generate/`.

The frozen portfolio is compatible with the AMP Challenge 2027 external evaluation protocol, including:

- 50,000 unique canonical peptide sequences;
- sequence lengths between 8 and 50 residues;
- no exact full-library overlap with the supplied antibacterial reference;
- a ranked 100-candidate subset contained within the 50,000-member library;
- an internal Top100 Levenshtein similarity safety gate of 0.78, stricter than the external protocol limit of 0.80;
- deterministic byte-identical output on repeated runs.

The external benchmark is one evaluation use case of the broader generation + pre-experimental screening workflow.
