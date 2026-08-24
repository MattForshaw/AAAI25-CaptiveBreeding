# Examples

Six small, synthetic scenarios. Each runs in under a second, reaches `OPTIMAL`, and
exercises one part of the input format so you can read a whole input file at a glance.

Identifiers are `E01`, `E02`, … and species are `Alpha`, `Beta` and `Gamma`. Relatedness comes from a fixed formula, `r(i,j) = 40 + 900 · ((7i + 13j) mod 11) / 10`.

| Example | Shows |
|---------|-------|
| [01-minimal](01-minimal/) | The required columns and a fixed-size group. Six individuals into three pairs — the smallest input that exercises the whole pipeline |
| [02-variable-size](02-variable-size/) | `MinSize` < `MaxSize`, a per-group `MaxPR` ceiling, and individuals left unallocated |
| [03-sex-ratios](03-sex-ratios/) | The optional `RatioMale`/`RatioFemale`/`RatioMode` columns, one group per mode |
| [04-multi-species](04-multi-species/) | `;`-separated species sets, an empty group meaning "any species", and `MaxNumNonComp` capping optional-species placements |
| [05-weighted-objective](05-weighted-objective/) | `WEIGHTED_ALLELES_PR_50_50` with component bounds derived from the data |
| [06-pre-allocation](06-pre-allocation/) | `AssignToFirstCorral` pinning an individual by zero-based group index, and `MustAllocate` |

Each directory has its own README with the exact command. To run them all:

```bash
for d in examples/0*/; do
  python3 cobreeder/solver.py run "$d/individuals.csv" "$d/pr-scaled.csv" \
      "$d/corrals.csv" ALL_PAIRS_PR_MIN_SQUARED "$(basename "$d")" 0 0 4 --quiet
done
```
