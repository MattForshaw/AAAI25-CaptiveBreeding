#!/usr/bin/env python3
"""Generate the example scenarios shipped in examples/.

These are **synthetic**. They exist so that someone can run the solver, exercise every
optional feature and read a small input file end to end without touching the study data.
Nothing here is from the paper: identifiers are E01, E02, ... and species are named Alpha,
Beta and Gamma precisely so no one can mistake them for the real cohorts.

Relatedness is generated from a fixed formula rather than a random draw, so the files are
byte-reproducible without carrying a seed:

    r(i, j) = base + spread * ((i * 7 + j * 13) mod 11) / 10

Committed so the examples can be regenerated and diffed:

    python3 examples/generate.py            # rewrite examples/*/
    python3 examples/generate.py --check    # verify the committed files match

Every example is small enough to solve in about a second, and every one is expected to
reach OPTIMAL. tests/test_examples.py asserts exactly that, so they cannot rot.
"""

from __future__ import annotations

import argparse
import csv
import io
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent

INDIVIDUAL_COLUMNS = ["Name", "Male", "Female", "AssignToFirstCorral", "Species",
                      "Alleles", "MustAllocate"]
CORRAL_COLUMNS = ["ID", "CompGroup", "OptionalGroup", "MinSize", "MaxSize",
                  "MinNumMale", "MinNumFemale", "MaxPR", "MaxNumNonComp"]


def relatedness(i, j, base=40, spread=900):
    """A deterministic stand-in for pairwise relatedness, on the solver's integer scale.

    Deliberately not random: the examples must be byte-reproducible, and a formula makes
    the values inspectable -- you can predict any cell by hand while reading the file.
    """
    if i == j:
        return 0
    a, b = min(i, j), max(i, j)
    return base + spread * ((a * 7 + b * 13) % 11) // 10


def render_csv(header, rows):
    """The CSV text for one file. Renders only -- main() decides whether to write it.

    It used to also create the destination directory, which meant --check, documented as
    "changing nothing", created every examples/*/ directory as a side effect.
    """
    buffer = io.StringIO()
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(header)
    writer.writerows(rows)
    return buffer.getvalue()


def individuals(n, *, sexes=None, species=None, alleles=None, must_allocate=None,
                first_corral=None):
    """One row per individual, in the same order as the relatedness matrix.

    Order matters: the solver pairs the two files by position, and validates that they
    agree, and rejects the file if they do not.
    """
    rows = []
    for i in range(n):
        male = (sexes[i] == "M") if sexes else (i % 2 == 1)
        rows.append([
            "E%02d (%s)" % (i + 1, "M" if male else "F"),
            1 if male else 0,
            0 if male else 1,
            (first_corral or {}).get(i, -1),
            (species or ["Alpha"] * n)[i],
            (alleles or [100] * n)[i],
            (must_allocate or [0] * n)[i],
        ])
    return INDIVIDUAL_COLUMNS, rows


def matrix(n, **kwargs):
    header = ["E%02d" % (i + 1) for i in range(n)]
    rows = [[relatedness(i, j, **kwargs) for j in range(n)] for i in range(n)]
    return header, rows


# --------------------------------------------------------------------- the examples


def example_minimal():
    """The smallest thing that solves: six individuals, three pairs."""
    n = 6
    files = {}
    files["individuals.csv"] = individuals(n)
    files["pr-scaled.csv"] = matrix(n)
    files["corrals.csv"] = (CORRAL_COLUMNS, [
        [i + 1, "Alpha", "", 2, 2, 1, 1, -1, -1] for i in range(3)
    ])
    return files, dict(
        title="Minimal",
        blurb="Six individuals into three mixed-sex pairs. The smallest input that "
              "exercises the whole pipeline.",
        objective="ALL_PAIRS_PR_MIN_SQUARED",
        total=6,
        shows="required columns, fixed-size groups, one male and one female per group",
    )


def example_variable_size():
    """Variable-size groups and a relatedness ceiling."""
    n = 12
    files = {}
    files["individuals.csv"] = individuals(n)
    files["pr-scaled.csv"] = matrix(n)
    # MaxPR bars any co-located pair above 500. MinSize < MaxSize is the case that was
    # the case a fixed-size-only configuration cannot express, so it is worth an example.
    files["corrals.csv"] = (CORRAL_COLUMNS, [
        [1, "Alpha", "", 3, 5, 1, 1, 500, -1],
        [2, "Alpha", "", 3, 5, 1, 1, 500, -1],
        [3, "Alpha", "", 2, 4, 1, 1, 500, -1],
    ])
    return files, dict(
        title="Variable-size groups",
        blurb="Groups that may hold a range of occupants, with a per-group ceiling on "
              "co-located relatedness.",
        objective="ALL_PAIRS_PR_MIN_SQUARED",
        total=10,
        shows="MinSize < MaxSize, MaxPR, and leaving individuals unallocated",
    )


def example_sex_ratio():
    """The optional ratio columns, one group per mode."""
    n = 12
    files = {}
    files["individuals.csv"] = individuals(
        n, sexes=["M", "M", "F", "F", "M", "F", "M", "F", "F", "M", "F", "F"])
    files["pr-scaled.csv"] = matrix(n)
    header = CORRAL_COLUMNS + ["RatioMale", "RatioFemale", "RatioMode"]
    files["corrals.csv"] = (header, [
        [1, "Alpha", "", 4, 4, 1, 1, -1, -1, 1, 1, "exact"],
        [2, "Alpha", "", 4, 4, 1, 1, -1, -1, 1, 2, "min-females"],
        [3, "Alpha", "", 4, 4, 1, 1, -1, -1, 1, 1, "min-males"],
    ])
    return files, dict(
        title="Sex ratios",
        blurb="One group per ratio mode: an exact 1:1, at least two females per male, "
              "and at least one male per female.",
        objective="ALL_PAIRS_PR_MIN_SQUARED",
        total=12,
        shows="the optional RatioMale/RatioFemale/RatioMode columns and all three modes",
    )


def example_multi_species():
    """Multi-species groups and a cap on optional-species placements."""
    n = 12
    species = (["Alpha"] * 6) + (["Beta"] * 4) + (["Gamma"] * 2)
    files = {}
    files["individuals.csv"] = individuals(
        n, species=species,
        sexes=["M", "F", "M", "F", "M", "F", "M", "F", "M", "F", "M", "F"])
    files["pr-scaled.csv"] = matrix(n)
    files["corrals.csv"] = (CORRAL_COLUMNS, [
        # Prefers Alpha, tolerates at most one Beta.
        [1, "Alpha", "Beta", 4, 4, 1, 1, -1, 1],
        # Treats Alpha and Beta equally, so neither counts as optional.
        [2, "Alpha;Beta", "", 4, 4, 1, 1, -1, -1],
        # Names no species, so it accepts anything -- including Gamma, which no other
        # group lists.
        [3, "", "", 4, 4, 1, 1, -1, -1],
    ])
    return files, dict(
        title="Multiple species",
        blurb="Three species across three groups: one preferring Alpha but tolerating a "
              "single Beta, one treating Alpha and Beta as equals, and one accepting "
              "anything.",
        objective="ALL_PAIRS_PR_MIN_SQUARED",
        total=12,
        shows="';'-separated species sets, an empty group meaning 'any species', and "
              "MaxNumNonComp capping optional-species placements",
    )


def example_weighted():
    """The two-component weighted objective."""
    n = 8
    files = {}
    # Varied allele counts, or the allele component has nothing to trade against.
    files["individuals.csv"] = individuals(
        n, alleles=[500, 120, 480, 90, 460, 140, 440, 110],
        must_allocate=[0] * n)
    files["pr-scaled.csv"] = matrix(n)
    files["corrals.csv"] = (CORRAL_COLUMNS, [
        [1, "Alpha", "", 2, 2, 1, 1, -1, -1],
        [2, "Alpha", "", 2, 2, 1, 1, -1, -1],
    ])
    return files, dict(
        title="Weighted objective",
        blurb="Trade total housed alleles against co-located relatedness. Vary the two "
              "weights to walk the trade-off; the ratio is what matters, not the "
              "magnitudes.",
        objective="WEIGHTED_ALLELES_PR_50_50",
        total=4,
        weights=(50, 50),
        shows="WEIGHTED_ALLELES_PR_50_50 with bounds derived from the data by "
              "--objective-bounds auto",
    )


def example_pre_allocation():
    """Pinning an individual to the first group, and requiring others be housed."""
    n = 9
    files = {}
    # first_corral values are zero-based *corral indices*, not flags: {0: 0} pins
    # individual 0 to corral 0. The column name reads like a yes/no about the first corral,
    # so 1 would pin the second group instead. tests/test_examples.py asserts the pinned
    # individual lands where the file names.
    files["individuals.csv"] = individuals(
        n, must_allocate=[1, 1, 1, 0, 0, 0, 0, 0, 0], first_corral={0: 0},
        sexes=["M", "F", "M", "F", "M", "F", "M", "F", "F"])
    files["pr-scaled.csv"] = matrix(n)
    files["corrals.csv"] = (CORRAL_COLUMNS, [
        [1, "Alpha", "", 3, 3, 1, 1, -1, -1],
        [2, "Alpha", "", 3, 3, 1, 1, -1, -1],
    ])
    return files, dict(
        title="Pre-allocation and required individuals",
        blurb="E01 is pinned to the first group; E01 to E03 must be housed; the rest are "
              "optional, and two of the nine are left out.",
        objective="ALL_PAIRS_PR_MIN_SQUARED",
        total=6,
        shows="AssignToFirstCorral pinning an individual to a group by zero-based index "
              "(despite the name it is not a flag), and MustAllocate=1 forcing placement "
              "while other individuals stay optional",
    )


EXAMPLES = {
    "01-minimal": example_minimal,
    "02-variable-size": example_variable_size,
    "03-sex-ratios": example_sex_ratio,
    "04-multi-species": example_multi_species,
    "05-weighted-objective": example_weighted,
    "06-pre-allocation": example_pre_allocation,
}


def command_for(name, meta):
    weights = meta.get("weights", (0, 0))
    return ("python3 cobreeder/solver.py run \\\n"
            "    examples/%s/individuals.csv \\\n"
            "    examples/%s/pr-scaled.csv \\\n"
            "    examples/%s/corrals.csv \\\n"
            "    %s %s %d %d %d" % (name, name, name, meta["objective"],
                                    name.split("-", 1)[1], weights[0], weights[1],
                                    meta["total"]))


def render(name, builder):
    files, meta = builder()
    out = {}
    for filename, (header, rows) in files.items():
        out[filename] = render_csv(header, rows)
    out["README.md"] = (
        "# %s\n\n%s\n\nSynthetic data. Identifiers are `E01`, `E02`, ... and species are "
        "`Alpha`, `Beta`, `Gamma`, so\nnothing here can be mistaken for the study "
        "cohorts.\n\n**Shows:** %s\n\n```bash\n%s\n```\n\nExpected outcome: `OPTIMAL`. "
        "Regenerate with `python3 examples/generate.py`.\n"
        % (meta["title"], meta["blurb"], meta["shows"], command_for(name, meta)))
    return out, meta


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--check", action="store_true",
                        help="verify the committed files match, changing nothing")
    args = parser.parse_args(argv)

    differences = []
    for name, builder in EXAMPLES.items():
        rendered, meta = render(name, builder)
        for filename, text in rendered.items():
            path = HERE / name / filename
            if args.check:
                current = path.read_text() if path.exists() else None
                if current != text:
                    differences.append(str(path.relative_to(REPO)))
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text(text)
        if not args.check:
            print("  examples/%s/  %s" % (name, meta["title"]))

    if args.check:
        if differences:
            print("out of date:", ", ".join(differences), file=sys.stderr)
            return 1
        print("all example files match examples/generate.py")
    return 0


if __name__ == "__main__":
    sys.exit(main())
