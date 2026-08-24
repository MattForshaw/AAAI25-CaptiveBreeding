# Pre-allocation and required individuals

E01 is pinned to the first group; E01 to E03 must be housed; the rest are optional, and two of the nine are left out.

Synthetic data. Identifiers are `E01`, `E02`, ... and species are `Alpha`, `Beta`, `Gamma`, so
nothing here can be mistaken for the study cohorts.

**Shows:** AssignToFirstCorral pinning an individual to a group by zero-based index (despite the name it is not a flag), and MustAllocate=1 forcing placement while other individuals stay optional

```bash
python3 cobreeder/solver.py run \
    examples/06-pre-allocation/individuals.csv \
    examples/06-pre-allocation/pr-scaled.csv \
    examples/06-pre-allocation/corrals.csv \
    ALL_PAIRS_PR_MIN_SQUARED pre-allocation 0 0 6
```

Expected outcome: `OPTIMAL`. Regenerate with `python3 examples/generate.py`.
