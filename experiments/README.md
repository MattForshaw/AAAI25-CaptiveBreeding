# Experiments

A parameterised runner over a committed manifest. One row per configuration, one command to
run any subset of them, and a second command to collect the results into CSVs.

```shell
experiments/run.sh --dry-run          # show what would run
experiments/run.sh a10-leaveout       # run one family
experiments/run.sh                    # run everything (slow, see Runtime below)
experiments/extract.sh                # collect results into CSVs
```

Output lands in `results/<run_id>-r<repeat>.out`, with developer traces alongside in
`.log`. Existing output is skipped unless `--force` is given, so an interrupted sweep
resumes by re-running the same command.

## The manifest

`manifest.csv` holds one row per configuration:

| Column | Meaning |
|--------|---------|
| `family` | Grouping label; `run.sh` takes these as filters |
| `run_id` | Stem for the output filename and the tag in every output record |
| `individuals`, `pr`, `corrals` | Repository-relative input paths |
| `objective` | Objective function name |
| `weight_a`, `weight_b` | Only used by `WEIGHTED_ALLELES_PR_50_50` |
| `total_individuals` | Minimum number of individuals allocated overall |
| `repeats` | How many times to run this configuration, each with `--seed <repeat>` |

`total_individuals` is a floor, not a target. For fixed-size groups the per-group size
constraints already imply it; for `pop-45to51`, where groups hold five or six, it is the
parameter being swept.

Add a row to add a run. Lines beginning `#` are ignored, and CRLF endings are tolerated so
the file can be edited in a spreadsheet.

## Families

| Family | Configs | Repeats | What it varies | Paper |
|--------|---------|---------|----------------|-------|
| `a10-leaveout` | 11 | 3 | `MaxPR` threshold, 10 individuals into 2 groups | — |
| `a20-leaveout` | 11 | 3 | `MaxPR` threshold, 20 individuals into 4 groups | — |
| `scaling` | 5 | 3 | Cohort size at a fixed group shape (groups of five) | Table 1, tortoise rows |
| `pop-45to51` | 7 | 3 | Individuals housed, 45→51, nine groups of 2M+3F | Fig 4 (centre) |
| `wolves` | 2 | 3 | Fixed-size vs variable-size groups, canid case study | Table 1, canid row |
| `wolves-weighted` | 21 | 1 | Objective weight α = 0, 0.05, … 1 | Fig 4 (right) |
| `multispecies` | 1 | 3 | Multi-species groups and a `MaxNumNonComp` cap | — |

## Runtime

132 runs in total, and they are not evenly matched:

| Family | Cost |
|--------|------|
| `a10-leaveout`, `wolves`, `wolves-weighted`, `multispecies` | seconds each |
| `scaling` at 10 and 20 individuals | under ten seconds each |
| `a20-leaveout` at the tighter thresholds | tens of minutes each |
| `scaling` at 30 individuals and above, `pop-45to51` | no first solution within several minutes |

Run families selectively. `analysis/figures.qmd` degrades gracefully when a family has not
been run, so it is worth rendering from a partial sweep.

## Where the families come from

Two are reconstructions of experiments the paper describes but the shipped data could not
express, and they are marked as such rather than presented as recovered configuration:

- **`pop-45to51`** raises the number of individuals housed from 45 to 51 within one cohort
  (paper §7.4, Figure 4 centre). No shipped group definition admits more than 45 — the
  nine-group files all set `MinSize = MaxSize = 5` — so
  `data/scenario-51indvs/pop45to51/corrals-pop45to51.csv` was written for this release:
  nine groups, `MinSize 5`, `MaxSize 6`, 2M+3F, no `MaxPR`.
- **`wolves-weighted`** sweeps α over 21 values, which is what §7.5 specifies: "α = 0, …,
  1 … in steps of 0.05".

Of the rest, `a10-leaveout` is reconstructed from a driver script that is not part of this
release — same eleven group files, same objective, same three repeats. Its input paths
named a directory this repository does not contain, so `data/scenario-a-10/` is used
instead, being the only cohort carrying those group files.

`a20-leaveout`, `scaling` and `wolves` are **inferred** from the group definitions in
`data/`. The groupings are the obvious reading of the filenames, but which rows correspond
to which published figure has not been confirmed against the authors' own records.
`objective` is `ALL_PAIRS_PR_MIN_SQUARED` throughout, the only objective named in the
surviving run instructions.

## Known gaps

Worth reading before treating the manifest as the paper's experiment set.

- **There is no 1000-individual scaling run.** That cohort is not part of this release: its
  relatedness matrix held 512,024 negative values, which the scaling step should have mapped
  onto a non-negative range as it does for the canids, and it carried no group definition.
  The `scaling` family therefore stops at 51 individuals. Reinstating it needs the matrix
  re-derived and a group definition supplied.
- **The 51-individual `MaxPR` percentile sweep is not included.** Its sixteen group files
  encoded percentiles p39–p54 of the relatedness distribution, but two of them carried the
  same ceiling, so the sweep never tested p42 or p43 and tested p44 three times — against a
  paper reporting p40–p45 as quickly infeasible. `corrals-realworld-maxpr.csv`, the
  operative threshold at p46, remains and is what `scaling` uses. Reinstating the sweep
  means re-deriving the percentiles and re-running.
- **Two `a10-leaveout` group files break the sweep's pattern.** Every other file sets group
  size `4,4`; files 14 and 15 set `4,10`. And `corrals-a10leaveout-15.csv` carries
  `MaxPR=22`, not 15, duplicating file 22 and leaving threshold 15 untested. The sweep
  therefore covers `MaxPR ∈ {12,13,14,16,…,22}` with 22 twice, rather than 12–22.
- **`MaxPR=12` is genuinely infeasible** on `a10-leaveout`, in all three repeats. That is a
  result, not a failure: no allocation of 8 of the 10 individuals into two groups of four
  keeps every co-located pair at or below relatedness 12.
- **`data/scenario-a/` and `data/scenario-a-50/` hold matrices only** and are not runnable;
  neither carries a cohort file in the current schema. `scenario-a` is retained because it
  also holds `updated_breeders_NGSrelate.csv`, a long-form pairwise source file.
- **`multispecies` is a demonstration, not a paper experiment.** `data/scenario-a-10/` is
  the only genuinely multi-species cohort (6 `R`, 4 `F`), and
  `corrals-multispecies.csv` was written for this release so the species and
  `MaxNumNonComp` features have a shipped configuration. A cap below 2 is infeasible there,
  because all four `F` individuals carry `MustAllocate=1` and two groups cannot admit them
  under a tighter cap.

### Two data questions that affect results

Both are recorded rather than resolved, because settling them is a decision for the study
authors. Neither has been applied to anything in `data/`.

- **Species labels.** Every tortoise cohort's `Species` column disagrees with
  `data/scenario-51indvs/51indvs_pairwise_ancestry.csv`, the long-form source, which records
  37 `R` and 14 `F` where the columns read `R` throughout. `MustAllocate` was derived from
  species, so where a family houses everybody this cannot bite, and where it leaves
  individuals out it decides *which*. Re-deriving species from the source and re-solving
  changes the proven optimum from 666 to 648 on `a10-leaveout` (thresholds 16, 18, 20, 22)
  and from 1557 to 1886 on `a20-leaveout-22`, while `scaling`, `multispecies` and `wolves`
  are unchanged — 66 of the 132 runs here. See `tests/test_data_provenance.py`.
- **Pre-allocation.** `AssignToFirstCorral` is a zero-based group index, and four cohorts
  set it: `scenario-51indvs` (two individuals), `scenario-a-20`, `scenario-a-30` and
  `scenario-a-40`. Two of them set `1`, which pins the *second* group. These are hard
  constraints on the `scaling`, `pop-45to51` and `a20-leaveout` families. On
  `scaling-scenario-a-20` the pin leaves the optimum at 6774 but changes the search from
  74,089 conflicts and 9 improving solutions to 52,284 and 14 — the quantities Table 1 and
  Figures 2 and 3 report.

## Reproducibility

The solver defaults to a single worker with a fixed seed, so a run repeats exactly.
`run.sh` passes `--seed <repeat number>`, which is what keeps the `repeats` column
meaningful: without it the repeats of a configuration would be byte-identical and measure
nothing. Each repeat is a distinct but individually reproducible sample of the search, so a
median across repeats is a meaningful summary.

`run.sh` writes `results/environment.txt` recording host, core count, Python and OR-Tools
versions, and the commit it ran from. Keep it with any results you keep.

Two caveats. A single worker is slower on hard instances than one per core — pass
`--workers 0` to the solver if you need speed more than reproducibility, accepting that the
result will differ between runs. And the published results were produced on a 64-vCPU
`c6i.16xlarge` with one worker per core, so their solve times are not comparable with
anything produced here.

## Extraction

`extract.sh` writes six CSVs into `results/`, each with a header row naming its columns:
`args.csv`, `solutions.csv`, `allocations.csv`, `completions.csv`, `search.csv` and
`bounds.csv`. The last two record what pinned each run — seed, workers, limits — and the
component ranges a weighted objective normalised against; `bounds.csv` is empty unless
`WEIGHTED_ALLELES_PR_50_50` was run.

Filenames are suppressed because every record already carries its `run_id`, so the files
load directly:

```r
completions <- read_csv("results/completions.csv")
```

`extract.sh` also checks every row's field count against its header and reports a mismatch
by name, which catches `.out` files produced by different solver versions being read
together.

### Reading the output with `AAAICR-Thesis.R`

That script supplies its own column names and expects the record tag in the first column,
which the default output does not provide: it strips the tag and adds a header. Use
`--legacy` for that consumer:

```shell
experiments/extract.sh --legacy
```

This keeps the tag as field 1 and writes no header, giving 10 columns for completions, 9
for solutions and 5 for allocations.
