# Multiple species

Three species across three groups: one preferring Alpha but tolerating a single Beta, one treating Alpha and Beta as equals, and one accepting anything.

Synthetic data. Identifiers are `E01`, `E02`, ... and species are `Alpha`, `Beta`, `Gamma`, so
nothing here can be mistaken for the study cohorts.

**Shows:** ';'-separated species sets, an empty group meaning 'any species', and MaxNumNonComp capping optional-species placements

```bash
python3 cobreeder/solver.py run \
    examples/04-multi-species/individuals.csv \
    examples/04-multi-species/pr-scaled.csv \
    examples/04-multi-species/corrals.csv \
    ALL_PAIRS_PR_MIN_SQUARED multi-species 0 0 12
```

Expected outcome: `OPTIMAL`. Regenerate with `python3 examples/generate.py`.
