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
