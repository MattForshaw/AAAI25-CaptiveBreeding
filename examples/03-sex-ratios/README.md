# Sex ratios

One group per ratio mode: an exact 1:1, at least two females per male, and at least one male per female.

Synthetic data. Identifiers are `E01`, `E02`, ... and species are `Alpha`, `Beta`, `Gamma`, so
nothing here can be mistaken for the study cohorts.

**Shows:** the optional RatioMale/RatioFemale/RatioMode columns and all three modes

```bash
python3 cobreeder/solver.py run \
    examples/03-sex-ratios/individuals.csv \
    examples/03-sex-ratios/pr-scaled.csv \
    examples/03-sex-ratios/corrals.csv \
    ALL_PAIRS_PR_MIN_SQUARED sex-ratios 0 0 12
```

Expected outcome: `OPTIMAL`. Regenerate with `python3 examples/generate.py`.
