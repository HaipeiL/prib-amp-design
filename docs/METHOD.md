# Method and ranking documentation

## 1. M2 generator

M2 is an adaptive edit Transformer with 86,536,746 parameters (12 Transformer layers, hidden dimension 768, 12 attention heads, FFN dimension 3072). The biological condition is target pMIC plus bacterial species. Content and cardinality conditioning are both enabled. Generation starts at length 8 and uses 20 iterative residue/edit steps with length bounds 8–50, delete scale 0.35, and insert scale 0.30.

The frozen portfolio run generated 249,856 raw candidates from q90 and q50 pMIC conditions for *E. coli*, *P. aeruginosa*, and *S. aureus*. No retraining was performed for the Challenge run.

## 2. R01 activity qualification

Candidates were required to use the 20 canonical amino acids and satisfy the external-evaluation length bounds. Duplicates and exact matches to the organizer antibacterial reference were removed. A frozen three-model independent pMIC ensemble scored the clean pool. R01 produced the 50,000-member reservoir and a provisional potency-oriented Top100.

## 3. PRIB experimental-readiness audit (R02)

The generator and evaluators were frozen before R02. A 12,000-candidate high-potential frontier was stress-tested with deterministic local counterfactuals:

- one alanine-like substitution at every position (`A`, or `G` when the original residue is `A`);
- one physicochemical-contrast substitution at every position;
- N-terminal truncation by 1 and 2 residues when the result remained at least 8 aa;
- C-terminal truncation by 1 and 2 residues under the same length constraint.

Every perturbation was re-scored by all three frozen evaluators across the three focal species.

R02 computed four percentile axes: baseline worst-species potency, perturbation lower-tail robustness, local retention within 0.5 pMIC, and evaluator consensus. The primary R02 reliability score is the minimum percentile across these axes. Sequence-only physicochemical descriptors are triage flags only; they are not toxicity predictions.

## 4. PRIB official-panel alignment audit (R04)

R04 does not generate new sequences and does not retrain any model. It re-scores the same frozen R02 12,000-candidate frontier with the same three evaluator checkpoints across the ten bacterial species represented in the organizer's 20-strain panel:

*Acinetobacter baumannii*, *Enterobacter cloacae*, *Escherichia coli*, *Klebsiella pneumoniae*, *Pseudomonas aeruginosa*, *Salmonella enterica*, *Bacillus subtilis*, *Staphylococcus aureus*, *Enterococcus faecalis*, and *Enterococcus faecium*.

Organizer strain multiplicities are used only as weights to summarize Broad, Gram-negative, Gram-positive, and MDR proxies. Because the evaluator conditions on species labels, **R04 is not a strain-specific MIC model** and does not distinguish individual ATCC/AIC/MDR isolates within the same species.

For each candidate, R04 computes percentile axes for:

1. broad-panel lower-tail predicted activity;
2. Gram-negative predicted activity;
3. Gram-positive predicted activity;
4. MDR-weighted predicted activity;
5. cross-model/species lower-tail activity;
6. evaluator consensus;
7. the frozen R02 reliability floor.

The final R04 reliability score is the **minimum** of these axes, followed by geometric mean as a secondary ranking term. This maximin rule prevents strong performance on one official category proxy from averaging away a severe weakness on another or a weak R02 robustness profile.

## 5. Reference and diversity gates

Candidates must retain R02 `PASS` status and have no severe physicochemical flag. High-ranking candidates are checked against the organizer antibacterial reference. The internal Top100 Levenshtein safety threshold is 0.78, stricter than the external-evaluation limit of 0.80. A 3-mer Jaccard diversity gate is then applied before filling the 100-member portfolio.

The final 50,000-member library is the final R04 Top100 followed by the frozen clean R02/R01 reservoir, preserving uniqueness and official sequence constraints.

## 6. Final selection result

R04 retained 50/100 of the R02 Top100 and replaced 50/100. All final 100 candidates remain R02 `PASS` and none carries a severe physicochemical flag.

Final median computational proxies:
- broad-panel lower tail: **5.3391**
- Gram-negative: **5.4957**
- Gram-positive: **5.5833**
- MDR: **5.5127**
- panel cross-model/species floor: **5.0673**
- retained R02 reliability floor: **0.6217**

These are computational proxies only. Experimental MIC, HC50, hemolysis, selectivity, and safety remain to be measured.
