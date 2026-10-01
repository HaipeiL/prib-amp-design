# PRIB Pan-Panel Audit R04

Status: **PASS**

R04 keeps the frozen M2 generator and R02 perturbation audit unchanged, then adds an independent species-level stress test aligned to the 10 unique bacterial species represented in the official 20-strain panel. The frozen three-model pMIC evaluator ensemble is queried for every candidate/species pair. Official strain multiplicities are used only as weights when summarizing Broad, Gram-negative, Gram-positive, and MDR activity proxies.

This is deliberately **not** presented as strain-specific prediction: the evaluator conditions on species labels, not ATCC/AIC/MDR strain identities. HC50/hemolysis is also not predicted. Physicochemical flags remain pre-experimental risk proxies only.

- R02 Top100 retained: 50/100
- R02 Top100 replaced: 50/100
- Median broad-panel lower-tail proxy: 5.3391
- Median Gram-negative proxy: 5.4957
- Median Gram-positive proxy: 5.5833
- Median MDR proxy: 5.5127
- Median cross-model/species floor: 5.0673
- Max official-reference similarity: 0.777778
