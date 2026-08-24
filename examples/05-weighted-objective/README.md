# Weighted objective

Trade total housed alleles against co-located relatedness. Vary the two weights to walk the trade-off; the ratio is what matters, not the magnitudes.

Synthetic data. Identifiers are `E01`, `E02`, ... and species are `Alpha`, `Beta`, `Gamma`, so
nothing here can be mistaken for the study cohorts.

**Shows:** WEIGHTED_ALLELES_PR_50_50 with bounds derived from the data by --objective-bounds auto

```bash
python3 cobreeder/solver.py run \
    examples/05-weighted-objective/individuals.csv \
    examples/05-weighted-objective/pr-scaled.csv \
    examples/05-weighted-objective/corrals.csv \
    WEIGHTED_ALLELES_PR_50_50 weighted-objective 50 50 4
```

Expected outcome: `OPTIMAL`. Regenerate with `python3 examples/generate.py`.
