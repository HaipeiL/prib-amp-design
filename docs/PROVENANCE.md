# Frozen provenance

## M2 generator checkpoint

- Stage: `finetune`
- Variant: `full`
- Step: `10000`
- Parameters: `86,536,746`
- Content conditioning: enabled
- Cardinality conditioning: enabled
- SHA256: `fdc17a1052dcd00678bbe317dbd7e0d29ee36ba36e5cc6221707f07ffb8e2141`

## Independent pMIC evaluator checkpoints

1. seed 20260903, step 8000 — SHA256 `a690836e8f9fe68e318a012fa35bb1c80615d608ea1173466ac6e49479eb181a`
2. seed 20260917, step 5000 — SHA256 `5eaf2024023c53fef9d5ba58ce3d375fc448a3eaad62cf3365ba6dd8f2ac1ea3`
3. seed 20261001, step 7000 — SHA256 `2881231b3a631ba874f0a9d6a287b864f9c99a2d882130a8238878cdf58e223b`

## Frozen portfolio generation

- Raw candidates: 249,856
- Unique/valid before exact-reference exclusion: 240,620
- Clean after exact-reference exclusion: 240,118
- Batch size: 128 per GPU
- Generation steps: 20
- Minimum/maximum length: 8/50
- Delete/insert scales: 0.35/0.30
- pMIC q50 raw: 5.008773924307505 (z = -0.012332823495805417)
- pMIC q90 raw: 5.946983986921648 (z = 1.154940058052101)

## PRIB audit

- Frontier: 12,000
- Eligible after audit + reference gate: 382
- Final Top100: 100
- R01 Top100 retained: 12
- R01 Top100 replaced: 88
- Maximum final Top100 reference similarity: 0.7777777777777778

## Source-of-truth rule

The final FASTA files in `assets/` are frozen portfolio artifacts. The audit CSV/JSON under `audit/` is the numerical source for the PRIB selection. `reference_code/` preserves the inference and audit implementation used during the Cluster42 run. Large model weights are tracked separately under the checkpoint export instructions below and must match the hashes above.


## R04 pan-panel finalization

- Input frontier: 12,000 frozen R02-audited candidates
- New training: none
- New generation: none
- Evaluators: same three frozen independent pMIC checkpoints listed above
- Species evaluated: 10 species underlying the official 20-strain panel
- Final Top100: 100
- R02 Top100 retained: 50
- R02 Top100 replaced: 50
- Final candidates retaining R02 PASS: 100/100
- Final candidates with severe physicochemical flag: 0/100
- Maximum final Top100 reference similarity: 0.7777777777777778
- Frozen final library SHA256: `0ece81746c2e514e870d68b38798c61669a56c6ed0df1b3e7084bc1d2a698a31`
- Frozen final Top100 SHA256: `623a3729bd57211466125d49bc90f89912ed275903b02e7489a020f66fd4e14a`

R04 is a species-conditioned computational proxy. It is not strain-specific prediction and does not model HC50/hemolysis.
