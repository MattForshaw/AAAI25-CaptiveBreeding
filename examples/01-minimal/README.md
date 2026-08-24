# Minimal

Six individuals into three mixed-sex pairs. The smallest input that exercises the whole pipeline.

Synthetic data. Identifiers are `E01`, `E02`, ... and species are `Alpha`, `Beta`, `Gamma`, so
nothing here can be mistaken for the study cohorts.

**Shows:** required columns, fixed-size groups, one male and one female per group

```bash
python3 cobreeder/solver.py run \
    examples/01-minimal/individuals.csv \
    examples/01-minimal/pr-scaled.csv \
    examples/01-minimal/corrals.csv \
    ALL_PAIRS_PR_MIN_SQUARED minimal 0 0 6
```

Expected outcome: `OPTIMAL`. Regenerate with `python3 examples/generate.py`.
