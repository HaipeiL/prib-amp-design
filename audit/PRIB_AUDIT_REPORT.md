# PRIB Experimental-Readiness Audit R02

Status: **PASS**

- Frontier stress-tested: 12000
- Eligible after PRIB + official reference gate: 382
- Final library: 50000
- Final top list: 100
- R01 Top100 retained: 12/100
- R01 Top100 replaced by the audit: 88/100
- Max reference similarity in final Top100: 0.777778
- Median broad-spectrum baseline pMIC proxy (min over 3 species): 5.3972
- Median perturbation lower-tail pMIC proxy (10th percentile): 4.9912
- Median cross-model perturbation lower-tail floor: 4.7605
- Median fraction retaining within 0.5 pMIC under perturbation: 0.929
- Median maximum evaluator disagreement SD: 0.1415

## Method
Generator output was frozen before this audit. R02 does not retrain M2. A 12k high-potential, six-stratum frontier is stress-tested using deterministic position-wise alanine-like substitutions, physicochemical-contrast substitutions, and N/C-terminal truncations. The same frozen independent pMIC ensemble is used only as a measurement instrument for local perturbation response. Selection is failure-oriented: potency, perturbation lower-tail performance, local retention, and evaluator consensus are converted to within-frontier percentiles; the candidate's audit score is the minimum (maximin) across those axes. Sequence-only physicochemical descriptors are used as caution/severe risk flags, not as toxicity claims. Official novelty/reference constraints are applied afterward, followed by a diversity gate.

## Interpretation boundary
This is a computational pre-experimental audit. It does not establish MIC, HC50, hemolysis, pharmacokinetics, or safety. Its purpose is to reduce model-fragile and physically extreme candidates before wet-lab selection.
