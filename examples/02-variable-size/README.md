# Variable-size groups

Groups that may hold a range of occupants, with a per-group ceiling on co-located relatedness.

Synthetic data. Identifiers are `E01`, `E02`, ... and species are `Alpha`, `Beta`, `Gamma`, so
nothing here can be mistaken for the study cohorts.

**Shows:** MinSize < MaxSize, MaxPR, and leaving individuals unallocated

```bash
python3 cobreeder/solver.py run \
    examples/02-variable-size/individuals.csv \
    examples/02-variable-size/pr-scaled.csv \
    examples/02-variable-size/corrals.csv \
    ALL_PAIRS_PR_MIN_SQUARED variable-size 0 0 10
```

Expected outcome: `OPTIMAL`. Regenerate with `python3 examples/generate.py`.
