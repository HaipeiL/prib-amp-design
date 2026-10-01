# Training data and external-data disclosure

## Generator pretraining

The common adaptive model was pretrained on a bank of **100,000 unique short canonical UniProt sequences**, spanning lengths 8–64 residues. The selected shared initialization was the 240,000-step pretraining checkpoint.

## Conditional AMP finetuning

M2 was finetuned on the **AMPBench-MT MIC** conditional task. Frozen project accounting records 50,660 training rows, 11,543 validation rows, and 11,453 test rows. The official AMPBench test split was **not used** for model selection, generation benchmarking, independent evaluator training, condition sweeps, or the frozen portfolio-selection workflow.

The biological condition is target pMIC plus bacterial species.

## Independent pMIC evaluators

Three frozen independent pMIC evaluator checkpoints were used as measurement instruments for R01 ranking and R02 perturbation response. They were trained within the frozen project protocol; the official AMPBench test split was not used for evaluator training or selection. The evaluators are not claimed to be state of the art and their outputs are treated as predicted pMIC proxies rather than experimental measurements.

## Organizer antibacterial reference

The organizer-provided `data/antibacterial.fasta` was used only as a post-generation reference gate:

- the 50,000-member library excludes exact sequence matches;
- the Top100 uses an internal maximum Levenshtein-ratio threshold of 0.78.

The organizer reference was not used to retrain M2 during the frozen portfolio run.

## Filters applied

- canonical 20-amino-acid alphabet;
- sequence length 8–50;
- uniqueness;
- exact organizer-reference exclusion for the library;
- independent predicted pMIC qualification across three focal species;
- PRIB perturbation robustness / evaluator-consensus audit;
- sequence-only physicochemical extreme-value flags;
- Top100 reference-similarity gate;
- Top100 3-mer diversity gate.
