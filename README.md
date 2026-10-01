# PRIB AMP Design

**Adaptive antimicrobial-peptide generation with pre-experimental robustness and panel-level screening.**

PRIB AMP Design is a research workflow for moving from computational peptide generation toward a more experimentally actionable candidate portfolio. It separates **candidate proposal** from **activity qualification**, **local perturbation stress testing**, and **multi-target panel screening**, rather than treating generator score alone as the final selection criterion.

## Workflow

```text
adaptive-cardinality generation
        ↓
independent activity qualification
        ↓
PRIB perturbation robustness audit
        ↓
multi-species panel audit
        ↓
novelty + diversity gates
        ↓
experimental candidate portfolio
```

### Adaptive-cardinality generation

M2 is an 86.54M-parameter adaptive edit Transformer conditioned on target pMIC and bacterial species. Residue content and sequence cardinality are both condition-dependent, allowing final peptide length to emerge through iterative insertion/deletion dynamics rather than being fixed in advance.

### Independent activity qualification

A frozen three-model pMIC ensemble scores the clean candidate pool independently of the generator.

### PRIB experimental-readiness audit

A 12,000-candidate frontier is subjected to deterministic local counterfactual stress tests, including position-wise substitutions and terminal truncations. Candidates are evaluated for lower-tail perturbation response, local activity retention, evaluator disagreement, and sequence-only physicochemical risk flags.

### Multi-species panel audit

The same frozen frontier is evaluated across ten clinically relevant bacterial species spanning Gram-negative and Gram-positive organisms. Selection uses a maximin criterion across broad-panel, subgroup, model-consensus, and retained perturbation-robustness axes so that a strong average score cannot hide a severe weakness.

## Frozen research portfolio

- Raw M2 generation: **249,856**
- Clean post-filter pool: **240,118**
- Perturbation-audited frontier: **12,000**
- Multi-species audited frontier: **12,000**
- Frozen design library: **50,000**
- Ranked experimental portfolio: **100**
- Final portfolio retaining PRIB audit PASS: **100/100**
- Severe physicochemical flags in final portfolio: **0/100**
- Internal maximum reference Levenshtein ratio: **0.777778**

These are **computational screening results**, not experimental MIC, strain-specific MIC, HC50, hemolysis, pharmacokinetic, or safety measurements.

## Deterministic reproduction

```bash
uv sync
uv run generate
```

This writes:

```text
generate/library.fasta
generate/top.fasta
```

The command materializes the frozen audited portfolio byte-for-byte so that an external benchmark or independent reviewer receives exactly the same candidate set.

Frozen SHA256:

```text
library: 0ece81746c2e514e870d68b38798c61669a56c6ed0df1b3e7084bc1d2a698a31
top100:  623a3729bd57211466125d49bc90f89912ed275903b02e7489a020f66fd4e14a
```

## Repository structure

```text
assets/          frozen 50k library and ranked 100-candidate portfolio
audit/           perturbation and multi-species audit evidence
checkpoint/      frozen model weights and hashes
docs/            method, data disclosure, provenance, evaluation notes
reference_code/  source used for generation, scoring, and audits
src/             deterministic reproduction interface
tests/           frozen-output integrity checks
```

## Documentation

- `docs/METHOD.md` — complete generator and screening workflow
- `docs/TRAINING_DATA.md` — training and external-data disclosure
- `docs/PROVENANCE.md` — checkpoint and artifact provenance
- `docs/EXTERNAL_EVALUATION.md` — deterministic external-benchmark interface
- `audit/` — complete perturbation and panel-audit evidence

## Scope

The PRIB layers are **pre-experimental computational audits**. They prioritize candidates for experimental validation; they do not establish efficacy, toxicity, selectivity, or safety.

## License

MIT.
