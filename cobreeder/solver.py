#!/usr/bin/env python3
import time
from typing import Sequence
from ortools.sat.python import cp_model
import pandas as pd
import argparse
from enum import IntEnum
import logging
import re
from fractions import Fraction
from math import gcd
import sys
from types import SimpleNamespace


INDIVIDUAL_COLUMNS = ["Name", "Male", "Female", "AssignToFirstCorral", "Species",
                      "Alleles", "MustAllocate"]

CORRAL_COLUMNS = ["ID", "CompGroup", "OptionalGroup", "MinSize", "MaxSize",
                  "MinNumMale", "MinNumFemale", "MaxPR", "MaxNumNonComp"]

# Optional columns. Absent from a corral file means the constraint is simply not
# applied, so existing files stay valid. RatioMale/RatioFemale must appear together;
# RatioMode is optional alongside them and defaults to DEFAULT_SEX_RATIO_MODE.
SEX_RATIO_COLUMNS = ["RatioMale", "RatioFemale"]
SEX_RATIO_MODE_COLUMN = "RatioMode"
SEX_RATIO_MODES = ("exact", "min-females", "min-males")
DEFAULT_SEX_RATIO_MODE = "exact"

# CP-SAT status -> (name for human-readable output, process exit code).
#
# A definitive answer exits 0 even when that answer is INFEASIBLE: in a threshold
# sweep, infeasibility is a result about the data, not a malfunction. The two
# non-zero codes distinguish the cases where the run produced no usable answer:
# MODEL_INVALID means the model was rejected before solving and is a defect, so it
# must stop a batch rather than quietly add a zero-solution row to a results file.
#
# Exit codes 1 and 2 are already in use -- 1 for input errors (UsageError) and
# uncaught exceptions, 2 by argparse for command-line misuse -- so this table
# starts at 3.
SOLVER_STATUS = {
    cp_model.OPTIMAL: ("OPTIMAL", 0),
    cp_model.FEASIBLE: ("FEASIBLE", 0),
    cp_model.INFEASIBLE: ("INFEASIBLE", 0),
    cp_model.MODEL_INVALID: ("MODEL_INVALID", 3),
    cp_model.UNKNOWN: ("UNKNOWN", 4),
}


INT64_MAX = 2 ** 63 - 1

# Developer traces go here, on stderr, so stdout carries only the run's results: the
# ###COBREEDER-* records that the analysis pipeline parses, the human-readable group
# listings, and the closing statistics. Before this split a 14-individual run put 738
# lines on stdout of which 152 were records.
__version__ = "0.1.0"          # keep in step with pyproject.toml

LOG = logging.getLogger("cobreeder")


def configure_logging(verbose=False, quiet=False):
    """Send traces to stderr. INFO by default, DEBUG with --verbose.

    DEBUG includes the per-cell corral/individual compatibility dump, the raw linear
    expressions, and CP-SAT's own search log -- together the bulk of the old output,
    and quadratic in population size.
    """
    level = logging.DEBUG if verbose else logging.WARNING if quiet else logging.INFO
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(logging.Formatter("%(levelname)s %(message)s"))
    LOG.handlers = [handler]
    LOG.setLevel(level)
    LOG.propagate = False


class UsageError(Exception):
    """A problem with the inputs or the command line.

    Reported as a message and exit status 1, without a traceback.
    """

# Components combined by the multi-objective mode, with the direction each is
# optimised in: +1 maximise, -1 minimise. Normalisation below is always
# increasing, and direction is applied only here, so the two cannot be confused.
WEIGHTED_COMPONENTS = (("alleles", +1), ("all_pairs_pr_squared", -1))


def reduce_sex_ratio(ratio_male, ratio_female, row, corral_file):
    """Reduce a male:female ratio to lowest terms, rejecting nonsensical values.

    Returns None when the ratio is disabled for this corral, otherwise (male, female)
    in lowest terms. Reducing matters because the ratio's period -- male + female --
    determines which occupancies can satisfy it, so 2:4 and 1:2 must behave alike.
    """
    if ratio_male < 0 or ratio_female < 0:
        return None
    if ratio_male == 0 and ratio_female == 0:
        raise UsageError(
            "Corral file %s row %d sets RatioMale and RatioFemale both to 0, which "
            "describes an empty corral. Use -1 in either column to disable the ratio."
            % (corral_file, row))
    divisor = gcd(ratio_male, ratio_female)
    return ratio_male // divisor, ratio_female // divisor


def sex_ratio_satisfied(ratio_male, ratio_female, mode, male_count, female_count):
    """Whether a (male, female) pair meets a ratio under the given mode.

    All three modes compare the same two products, so each is one linear constraint
    differing only in its operator:

      exact        f*males == m*females   the corral sits exactly on the ratio
      min-females  f*males <= m*females   at least f females for every m males
      min-males    f*males >= m*females   at least m males for every f females

    'exact' is the strictest reading of "a particular ratio". The two inequality
    modes express the commoner husbandry requirement -- enough females per male, say
    -- and, unlike 'exact', place no restriction on corral size.
    """
    left = ratio_female * male_count
    right = ratio_male * female_count
    if mode == "exact":
        return left == right
    if mode == "min-females":
        return left <= right
    if mode == "min-males":
        return left >= right
    raise UsageError("Unknown sex-ratio mode %r" % mode)


def sex_ratio_occupancies(ratio_male, ratio_female, mode, min_size, max_size,
                          num_male, num_female):
    """Occupancies that can satisfy a ratio for one corral, under its mode.

    Enumerated rather than derived because the answer depends on the mode. Under
    'exact' only whole multiples of m + f qualify, since the corral holds k*m males
    and k*f females; the inequality modes are far looser. Either way the occupancy
    must also clear the MinNumMale and MinNumFemale floors, which is where a ratio and
    those minimums can silently disagree.
    """
    # Clamped at zero for the same reason min_size is: a negative floor means "no floor",
    # and an unclamped range would enumerate a negative count of one sex.
    num_male, num_female = max(num_male, 0), max(num_female, 0)
    feasible = []
    for occupancy in range(max(min_size, 0), max_size + 1):
        for male_count in range(num_male, occupancy - num_female + 1):
            female_count = occupancy - male_count
            if sex_ratio_satisfied(ratio_male, ratio_female, mode,
                                   male_count, female_count):
                feasible.append(occupancy)
                break
    return feasible


def validate_sex_ratios(corral_defs, corral_file):
    """Reject ratios no allocation could satisfy, naming the row and the reason.

    Without this an unsatisfiable ratio surfaces only as INFEASIBLE for the whole
    model, which says nothing about which corral is at fault.

    Returns {corral_index: (male, female, mode)} for the corrals with a ratio set.
    """
    if not all(c in corral_defs.columns for c in SEX_RATIO_COLUMNS):
        missing = [c for c in SEX_RATIO_COLUMNS if c not in corral_defs.columns]
        present = [c for c in SEX_RATIO_COLUMNS if c in corral_defs.columns]
        if present:
            raise UsageError(
                "Corral file %s has %s but not %s. The sex-ratio columns must be "
                "supplied together." % (corral_file, ", ".join(present),
                                        ", ".join(missing)))
        return {}

    has_mode = SEX_RATIO_MODE_COLUMN in corral_defs.columns
    ratios = {}
    for t in range(len(corral_defs)):
        reduced = reduce_sex_ratio(int(corral_defs["RatioMale"][t]),
                                   int(corral_defs["RatioFemale"][t]),
                                   t, corral_file)
        if reduced is None:
            continue
        ratio_male, ratio_female = reduced

        mode = DEFAULT_SEX_RATIO_MODE
        if has_mode:
            raw = corral_defs[SEX_RATIO_MODE_COLUMN][t]
            mode = DEFAULT_SEX_RATIO_MODE if pd.isna(raw) else str(raw).strip()
            if mode not in SEX_RATIO_MODES:
                raise UsageError(
                    "Corral file %s row %d sets RatioMode to %r. Valid modes are %s."
                    % (corral_file, t, mode, ", ".join(SEX_RATIO_MODES)))

        min_size = int(corral_defs["MinSize"][t])
        max_size = int(corral_defs["MaxSize"][t])
        num_male = int(corral_defs["MinNumMale"][t])
        num_female = int(corral_defs["MinNumFemale"][t])

        feasible = sex_ratio_occupancies(ratio_male, ratio_female, mode,
                                         min_size, max_size, num_male, num_female)
        if not feasible:
            detail = ""
            if mode == "exact":
                period = ratio_male + ratio_female
                # Show a few multiples even when none reach max_size, so the message
                # says what the ratio does allow rather than only what it does not.
                multiples = [k * period
                             for k in range(1, max(4, max_size // period + 2))]
                detail = (" An exact ratio only admits occupancies that are multiples "
                          "of %d (%s)." % (period,
                                           ", ".join(str(m) for m in multiples[:6])))
            raise UsageError(
                "Corral file %s row %d asks for a %d:%d male:female ratio in %s mode, "
                "which no allocation can satisfy at size %d-%d while also meeting "
                "MinNumMale >= %d and MinNumFemale >= %d.%s Widen the size range, relax the "
                "minimums, choose another ratio, or use a less strict RatioMode."
                % (corral_file, t, ratio_male, ratio_female, mode,
                   min_size, max_size, num_male, num_female, detail))
        ratios[t] = (ratio_male, ratio_female, mode)
    return ratios


def individual_identifiers(names):
    """The identifier portion of each individual's Name.

    Name conventions differ between scenarios -- "LA34M CL14229 (Cameron, Moore-Odom,
    Transient)", "FloHyb_26055848 (Sex=F, Species=R)", "FloHyb_26055848 (F) (Spec:F)" --
    but every one leads with the identifier, and that identifier is what the relatedness
    matrix uses as its column heading. So the comparison is on the leading
    whitespace-delimited token, not on the whole name.

    Token equality rather than a prefix test, deliberately: several cohorts contain
    identifiers that are proper prefixes of others, so startswith() would accept a
    misalignment. The same rule lives in preprocessing/scale_relatedness.py.
    """
    return [str(name).strip().strip('"').split(" ")[0] for name in names]


def validate_relatedness(pr, names, path):
    """Check the relatedness matrix against the cohort it is supposed to describe.

    The matrix is consumed positionally -- connections[g1][g2] -- so a matrix of the right
    size whose rows are in a different order from individuals.csv optimises the wrong
    values and returns a plausible, wrong answer with no error anywhere. Size alone was
    the only thing previously checked.

    Returns the matrix as a list of rows.
    """
    identifiers = individual_identifiers(names)

    if len(pr.index) != len(names):
        raise UsageError(
            "Relatedness file %s describes %d individuals but individuals file lists %d. "
            "The matrix is read positionally, one row per individual in file order."
            % (path, len(pr.index), len(names)))
    if len(pr.columns) != len(names):
        raise UsageError(
            "Relatedness file %s has %d columns but %d rows are expected; the matrix must "
            "be square, with one header cell per individual."
            % (path, len(pr.columns), len(names)))

    # Ordering. pandas suffixes duplicate headers (.1, .2), so strip that before
    # comparing -- a duplicated identifier is reported below rather than mistaken for a
    # mismatch.
    headers = [re.sub(r"\.\d+$", "", str(c)).strip() for c in pr.columns]
    shared = set(headers) & set(identifiers)
    if shared:
        if headers != identifiers:
            if sorted(headers) == sorted(identifiers):
                first = next(i for i, (h, n) in enumerate(zip(headers, identifiers))
                             if h != n)
                raise UsageError(
                    "Relatedness file %s lists the same individuals as the individuals "
                    "file but in a different order: row %d is %r in the matrix and %r in "
                    "individuals.csv. The matrix is read positionally, so this would "
                    "optimise the wrong relatedness values and still report a solution. "
                    "Reorder one file to match the other."
                    % (path, first, headers[first], identifiers[first]))
            only_matrix = [h for h in headers if h not in set(identifiers)]
            only_cohort = [n for n in identifiers if n not in set(headers)]
            raise UsageError(
                "Relatedness file %s does not describe the same individuals as the "
                "individuals file. In the matrix only: %s. In individuals.csv only: %s."
                % (path, ", ".join(only_matrix[:5]) or "none",
                   ", ".join(only_cohort[:5]) or "none"))
    else:
        # A header carrying no identifiers is unusual but not wrong -- a matrix written
        # without names gets V1..Vn or 0..n-1. Ordering cannot be checked in that case,
        # so say so instead of guessing.
        LOG.warning(
            "Relatedness file %s has no column headings matching the individuals file, so "
            "row ordering could not be verified. The matrix is read positionally: row i "
            "must be the individual on line i+1 of %s.", path, "individuals.csv")

    values = pr.values

    # Values must be non-negative: the objectives sum coefficient * colocated over pairs,
    # so a negative coefficient would reward co-locating that pair, and the old code
    # instead skipped such pairs entirely -- silently excluding them from the objective.
    # Scaling is expected to map raw relatedness onto a non-negative range before it
    # reaches the solver.
    negative = (values < 0)
    if negative.any():
        raise UsageError(
            "Relatedness file %s holds %d negative values, the lowest %s. Relatedness "
            "reaching the solver must be zero or positive; rescale the matrix onto a "
            "non-negative range first. Note a min-max rescale preserves the ordering "
            "that clamping to zero would destroy, so a least-related pair still scores "
            "below an unrelated one. preprocessing/scale_relatedness.py --mode minmax "
            "does this." % (path, int(negative.sum()), values.min()))

    # Symmetry is required, not cosmetic: every objective and constraint reads
    # connections[g1][g2] with g1 < g2, so the lower triangle is never consulted. An
    # asymmetric matrix therefore has half its data silently ignored, and the two halves
    # disagreeing about a pair means one of them is wrong.
    asymmetry = abs(values - values.T)
    if asymmetry.max() > 1e-9:
        worst = divmod(int(asymmetry.argmax()), asymmetry.shape[1])
        raise UsageError(
            "Relatedness file %s is not symmetric: [%d][%d] is %s but [%d][%d] is %s. "
            "Only the upper triangle is read, so the lower half would be ignored rather "
            "than averaged."
            % (path, worst[0], worst[1], values[worst[0]][worst[1]],
               worst[1], worst[0], values[worst[1]][worst[0]]))

    # The diagonal is never read either, so a non-zero one cannot change a result -- and
    # self-relatedness of 1 is a legitimate convention in some kinship matrices. Warn
    # rather than reject, since it usually means self-pairs leaked into the matrix build.
    diagonal = values.diagonal()
    if (abs(diagonal) > 1e-9).any():
        LOG.warning(
            "Relatedness file %s has %d non-zero value(s) on its diagonal (largest %s), "
            "which would make an individual related to itself. The diagonal is never "
            "read, so results are unaffected, but check how the matrix was built.",
            path, int((abs(diagonal) > 1e-9).sum()), diagonal.max())

    return values.tolist()


# Several cohorts embed the species in the Name as well as carrying it in the Species
# column -- "FloHyb_26056533 (Sex=M, Species=F)" and "FloHyb_26055848 (F) (Spec:F)". Only
# the column reaches the model, so a disagreement between the two is invisible at run
# time while changing which corrals an individual may enter and, where MustAllocate was
# derived from species, who must be housed at all.
SPECIES_ANNOTATION = re.compile(r"Spec(?:ies)?[:=]\s*([^,)\s]+)")


def check_species_annotations(names, species, path):
    """Warn where a Name's embedded species annotation contradicts the Species column.

    Advisory, not fatal: the annotation is free text and the column is what the model
    reads, so rejecting the file would refuse data that solves. But the two disagreeing
    means one of them is wrong, and in this repository both cases turned out to be the
    column -- see the species-label note in the README.

    Returns the disagreements as (index, identifier, annotated, column) tuples.
    """
    disagreements = []
    for g, name in enumerate(names):
        match = SPECIES_ANNOTATION.search(str(name))
        if match and match.group(1) != str(species[g]):
            disagreements.append((g, individual_identifiers([name])[0],
                                  match.group(1), str(species[g])))
    if disagreements:
        shown = ", ".join("row %d %s (name says %r, column says %r)" % d
                          for d in disagreements[:5])
        LOG.warning(
            "Individuals file %s: %d of %d rows carry a species in Name that "
            "contradicts the Species column: %s%s. Only the column reaches the model. "
            "Check both against the source genetics before trusting a result.",
            path, len(disagreements), len(names), shown,
            ", ..." if len(disagreements) > 5 else "")
    return disagreements


# CompGroup and OptionalGroup each hold zero or more species codes. The file is CSV, so
# the separator cannot be a comma; ';' is canonical and '|' accepted, because a corral
# file edited in a spreadsheet tends to acquire whichever the author reached for first.
SPECIES_SEPARATORS = ";|"


def parse_species_set(cell):
    """Zero or more species codes from one CompGroup/OptionalGroup cell.

    An empty cell -- which pandas reads as NaN, as data/wolves/corrals-a.csv does for
    OptionalGroup -- gives the empty set. Codes are compared case-sensitively, matching
    how the Species column in individuals.csv has always been read.
    """
    if cell is None:
        return frozenset()
    text = str(cell).strip()
    if not text or text.lower() == "nan":
        return frozenset()
    for sep in SPECIES_SEPARATORS[1:]:
        text = text.replace(sep, SPECIES_SEPARATORS[0])
    return frozenset(part.strip() for part in text.split(SPECIES_SEPARATORS[0])
                     if part.strip())


def validate_species_groups(corral_defs, corral_file, species, must_allocate,
                            individuals_file):
    """Resolve each corral's preferred and optional species, rejecting what cannot work.

    Returns {corral_index: (preferred, optional)} as frozensets, with `optional` already
    stripped of anything also preferred -- a species named in both is preferred, so it
    must not also count towards the MaxNumNonComp cap.

    Semantics, which were previously only implicit and unenforced:

      * An individual may be placed in a corral only if its species is preferred or
        optional there.
      * A corral naming no species at all accepts any species. That is the reading that
        keeps existing files working, and the alternative -- no species means nobody
        qualifies -- would make such a corral trivially infeasible.
      * MaxNumNonComp caps how many of a corral's occupants come from the optional set
        rather than the preferred one. -1 means no cap.
    """
    cohort = set(species)
    groups = {}
    unknown_declared = {}

    for t in range(len(corral_defs)):
        preferred = parse_species_set(corral_defs["CompGroup"][t])
        optional = parse_species_set(corral_defs["OptionalGroup"][t])

        both = preferred & optional
        if both:
            LOG.warning(
                "Corral file %s row %d names %s in both CompGroup and OptionalGroup. "
                "Treating it as preferred, so it does not count towards MaxNumNonComp.",
                corral_file, t, ", ".join(sorted(both)))
        optional = optional - preferred

        for declared in (preferred | optional) - cohort:
            unknown_declared.setdefault(declared, []).append(t)

        cap = corral_defs["MaxNumNonComp"][t]
        if cap != -1:
            if cap < 0:
                raise UsageError(
                    "Corral file %s row %d has MaxNumNonComp %s. Use -1 for no cap, or "
                    "zero or more to cap optional-species placements."
                    % (corral_file, t, cap))
            if not optional:
                LOG.warning(
                    "Corral file %s row %d caps optional-species placements at %d but "
                    "declares no optional species beyond its preferred ones, so the cap "
                    "cannot bind.", corral_file, t, int(cap))

        groups[t] = (preferred, optional)

    # A species declared by a corral but held by nobody is harmless and shipped:
    # every scenario-51indvs corral declares OptionalGroup=F against an all-R cohort.
    # Say so once rather than per row, because it usually means a copied corral file.
    for declared, rows in sorted(unknown_declared.items()):
        LOG.warning(
            "Corral file %s declares species %r (rows %s) that no individual in %s "
            "has, so it can never be matched.", corral_file, declared,
            ", ".join(str(r) for r in rows[:6]) + ("..." if len(rows) > 6 else ""),
            individuals_file)

    # An individual that must be housed but is welcome nowhere makes the whole model
    # infeasible, which on its own says nothing about why. Name it here instead.
    accepts_anything = any(not pref and not opt for pref, opt in groups.values())
    if not accepts_anything:
        welcome = set()
        for preferred, optional in groups.values():
            welcome |= preferred | optional
        stranded = sorted({species[g] for g in range(len(species))
                           if species[g] not in welcome and int(must_allocate[g]) == 1})
        if stranded:
            raise UsageError(
                "Species %s in %s must be allocated (MustAllocate=1) but no corral in %s "
                "lists %s as preferred or optional, so no allocation exists. Add the "
                "species to a corral's CompGroup or OptionalGroup, or set MustAllocate=0."
                % (", ".join(repr(x) for x in stranded), individuals_file, corral_file,
                   "it" if len(stranded) == 1 else "them"))
        optional_only = sorted({species[g] for g in range(len(species))
                                if species[g] not in welcome})
        if optional_only:
            LOG.warning(
                "Species %s in %s is listed by no corral in %s, so those individuals "
                "cannot be placed anywhere. They are left unallocated.",
                ", ".join(repr(x) for x in optional_only), individuals_file, corral_file)

    return groups


def sex_ratio_placement_cap(sex_ratios, num_corrals, cohort_males, cohort_females):
    """How many individuals a uniform sex ratio permits to be placed at all.

    A ratio of m:f in every corral forces f*males == m*females across the allocated
    set, so the scarcer sex caps the total: at most (m+f) * min(males//m, females//f)
    individuals can be housed however many corrals there are. With a 6-male,
    8-female cohort and a 1:1 ratio that is 12, not 14 -- a limit that otherwise
    shows up only as a bare INFEASIBLE.

    Returns None unless every corral carries the same ratio in 'exact' mode, since
    mixed ratios have no single closed form worth guessing at and the inequality modes
    impose no such cap.
    """
    if len(sex_ratios) != num_corrals or not sex_ratios:
        return None
    distinct = set(sex_ratios.values())
    if len(distinct) != 1:
        return None
    ratio_male, ratio_female, mode = distinct.pop()
    if mode != "exact":
        # Only an equality pins the aggregate composition. The inequality modes leave
        # the scarcer sex free to be topped up by the other, so no such cap exists.
        return None
    groups = []
    if ratio_male:
        groups.append(cohort_males // ratio_male)
    if ratio_female:
        groups.append(cohort_females // ratio_female)
    if not groups:
        return None
    return (ratio_male + ratio_female) * min(groups)


def normalise_terms(terms, lo, hi, weight, scale, name):
    """Min-max normalise a linear component, exactly, in integer arithmetic.

    A component is sum(coeff * bool_var). CP-SAT cannot divide an expression by a
    constant, so the division that min-max normalisation requires is applied to the
    coefficients instead -- they are ordinary Python integers at build time.

    Scaling each coefficient by weight*scale/(hi-lo) makes the component span about
    [0, weight*scale] as its raw value moves from lo to hi, so two components with
    wildly different magnitudes become directly comparable and `scale` bounds the
    objective's magnitude regardless of the data.

    Coefficients round individually, so the achievable span drifts from weight*scale
    by at most half the number of terms. That sets the precision floor on `scale`:
    1e6 leaves drift under 0.1% for the 1275 pairs of a 51-individual scenario, but a
    population large enough to give hundreds of thousands of pairs needs a
    proportionally larger scale -- roughly 200 times the pair count keeps drift under
    0.25%.

    Returns (scaled_terms, offset, magnitude_bound).
    """
    span = hi - lo
    if span <= 0:
        raise UsageError(
            "Component %s has bounds [%d, %d], giving a range of %d. Normalisation "
            "needs a positive range; supply --objective-bounds explicitly."
            % (name, lo, hi, span))

    factor = Fraction(weight * scale, span)
    scaled = [(int(round(coeff * factor)), var) for coeff, var in terms]
    scaled = [(coeff, var) for coeff, var in scaled if coeff != 0]
    offset = int(round(lo * factor))
    magnitude = sum(abs(coeff) for coeff, _ in scaled) + abs(offset)
    return scaled, offset, magnitude


def derive_component_bounds(model, solver_factory, components, time_limit):
    """Find each component's achievable range on the fully-constrained model.

    Deriving bounds rather than hardcoding them keeps scenario-specific values out
    of the source and makes the normalisation correct for whatever data is supplied.

    Both directions must prove optimality: a merely feasible answer is only a bound
    on the true optimum, and feeding that in as lo or hi silently corrupts the range
    -- which is the weighting.
    """
    bounds = {}
    for name, expression in components.items():
        for sense in ("min", "max"):
            model.ClearObjective()
            if sense == "min":
                model.Minimize(expression)
            else:
                model.Maximize(expression)
            solver = solver_factory()
            if time_limit:
                solver.parameters.max_time_in_seconds = time_limit
            status = solver.Solve(model)
            if status != cp_model.OPTIMAL:
                status_name = SOLVER_STATUS.get(status, ("status %i" % status, 0))[0]
                raise UsageError(
                    "Could not derive the %s bound for component %s: the solver "
                    "returned %s rather than OPTIMAL. A non-optimal bound would "
                    "distort the weighting, so supply --objective-bounds "
                    "LO_A,HI_A,LO_B,HI_B instead, or raise --bounds-time-limit."
                    % (sense, name, status_name))
            bounds[(name, sense)] = int(round(solver.ObjectiveValue()))
    model.ClearObjective()
    return bounds


class CobreederObjectiveFunction(IntEnum):
    ALL_PAIRS = 1
    MALE_FEMALE = 2
    ALLELES_MIN = 3
    ALLELES_MAX = 4
    ALL_PAIRS_PR_MIN = 5
    ALL_PAIRS_PR_MAX = 6
    WEIGHTED_ALLELES_PR_50_50 = 7
    ALL_PAIRS_PR_MIN_SQUARED = 8
    ALL_PAIRS_PR_MAX_SQUARED = 9
    MALE_FEMALE_SQUARED = 10


def terms_expression(terms):
    """Sum a list of (coefficient, variable) pairs into a linear expression."""
    return sum(coefficient * variable for coefficient, variable in terms)


def objective_alleles(ctx):
    """Total ghost-allele count over allocated individuals."""
    return terms_expression(ctx.allele_terms)


def objective_all_pairs_pr(ctx):
    """Summed relatedness over co-located pairs."""
    return sum(ctx.connections[a][b] * ctx.colocated[a, b] for a, b in ctx.pairs)


def objective_all_pairs_pr_squared(ctx):
    """Summed squared relatedness, which penalises close kin disproportionately."""
    return terms_expression(ctx.prsq_terms)


def objective_opposing_sex_pr(ctx):
    """Summed relatedness over opposite-sex co-located pairs only."""
    return sum(ctx.connections[a][b] * ctx.colocated[a, b]
               for a, b in ctx.pairs if ctx.opposing_sex[a, b])


def objective_opposing_sex_pr_squared(ctx):
    """Summed squared relatedness over opposite-sex co-located pairs only."""
    return sum(ctx.connections[a][b] * ctx.connections[a][b] * ctx.colocated[a, b]
               for a, b in ctx.pairs if ctx.opposing_sex[a, b])


# objective -> (sense, builder). Only the selected builder runs. Every pair-based
# expression is quadratic in population size, so constructing all of them on every run
# -- as the previous if/elif chain did -- wasted most of that work: on the
# 1000-individual scenario each is roughly half a million terms.
OBJECTIVE_BUILDERS = {
    CobreederObjectiveFunction.ALL_PAIRS: ("min", objective_all_pairs_pr),
    CobreederObjectiveFunction.ALL_PAIRS_PR_MIN: ("min", objective_all_pairs_pr),
    CobreederObjectiveFunction.ALL_PAIRS_PR_MAX: ("max", objective_all_pairs_pr),
    CobreederObjectiveFunction.ALL_PAIRS_PR_MIN_SQUARED: ("min", objective_all_pairs_pr_squared),
    CobreederObjectiveFunction.ALL_PAIRS_PR_MAX_SQUARED: ("max", objective_all_pairs_pr_squared),
    CobreederObjectiveFunction.MALE_FEMALE: ("min", objective_opposing_sex_pr),
    CobreederObjectiveFunction.MALE_FEMALE_SQUARED: ("min", objective_opposing_sex_pr_squared),
    CobreederObjectiveFunction.ALLELES_MIN: ("min", objective_alleles),
    CobreederObjectiveFunction.ALLELES_MAX: ("max", objective_alleles),
}

# Components the multi-objective mode normalises and combines, by name.
WEIGHTED_COMPONENT_BUILDERS = {
    "alleles": objective_alleles,
    "all_pairs_pr_squared": objective_all_pairs_pr_squared,
}


class CobreederPrinter(cp_model.CpSolverSolutionCallback):
    """Print intermediate solutions."""

    def __init__(self, placement, names, num_corrals, num_individuals, paramstring, unique_id):
        cp_model.CpSolverSolutionCallback.__init__(self)
        self.__solution_count = 0
        self.__start_time = time.time()
        self.__placement = placement
        self.__names = names
        self.__num_corrals = num_corrals
        self.__num_individuals = num_individuals
        self.__paramstring = paramstring
        self.__uniqueid = unique_id

    def on_solution_callback(self):
        current_time = time.time()
        objective = self.ObjectiveValue()

        self.__solution_count += 1

        print(
            "###COBREEDER-SOLUTION,%i,%f,%f,%i,%s,%s"
            % (self.__solution_count,
               current_time,
               current_time - self.__start_time,
               objective, self.__paramstring, self.__uniqueid)
        )

        for t in range(self.__num_corrals):
            print("Corral %d: " % t)
            for g in range(self.__num_individuals):
                if self.Value(self.__placement[(t, g)]):
                    print("  " + self.__names[g])

        for g in range(self.__num_individuals):
            print("###COBREEDER-ALLOCATION,%d,%d,%d,%s" % (self.__solution_count, g, self.get_corral_number(g), self.__uniqueid))

    def num_solutions(self):
        return self.__solution_count

    def get_corral_number(self, g):
        val = None
        for t in range(self.__num_corrals):
            if self.Value(self.__placement[(t, g)]):
                val = t
        return -1 if val is None else val


def build_data(args):
    """Build the data model."""
    objective_function = CobreederObjectiveFunction[args.obj_function]

    unique_id = args.unique_run_id

    def load(path, what):
        """Read an input file, reporting a missing or malformed one by name."""
        try:
            return pd.read_csv(path, delimiter=',')
        except FileNotFoundError:
            raise UsageError("No such %s file: %s" % (what, path))
        except IsADirectoryError:
            raise UsageError("%s is a directory, not a %s file" % (path, what))
        except PermissionError:
            raise UsageError("Cannot read %s file %s: permission denied" % (what, path))
        except pd.errors.EmptyDataError:
            raise UsageError("%s file %s is empty" % (what.capitalize(), path))
        except pd.errors.ParserError as exc:
            raise UsageError("%s file %s is not valid CSV: %s" % (what.capitalize(), path, exc))

    individuals = load(args.individuals_file, "individuals")

    # Checked before any column is read, so a missing one is explained rather than
    # surfacing as a bare pandas KeyError.
    missing_individual_columns = [c for c in INDIVIDUAL_COLUMNS
                                  if c not in individuals.columns]
    if missing_individual_columns:
        raise UsageError(
            "Individuals file %s is missing required column(s): %s. Files without "
            "Species and Alleles predate those columns. MustAllocate replaced an "
            "inference from Species: set it to 1 where an individual must be housed "
            "and 0 where placement is optional."
            % (args.individuals_file, ", ".join(missing_individual_columns)))

    names = individuals["Name"].tolist()
    males = individuals["Male"].tolist()
    females = individuals["Female"].tolist()
    allocate_first_corral = individuals["AssignToFirstCorral"].tolist()
    species = individuals["Species"].tolist()
    allele_counts = individuals["Alleles"].tolist()
    must_allocate = individuals["MustAllocate"].tolist()

    check_species_annotations(names, species, args.individuals_file)

    bad_must_allocate = sorted(set(int(v) for v in must_allocate) - {0, 1})
    if bad_must_allocate:
        raise UsageError(
            "Individuals file %s has MustAllocate values %s. Only 0 (placement "
            "optional) and 1 (must be housed) are meaningful."
            % (args.individuals_file, bad_must_allocate))

    pr = load(args.pairwise_relatedness_file, "relatedness")
    corral_defs = load(args.corral_file, "corral")

    # Reject unusable corral definitions here rather than failing part-way through
    # model building. Files carrying a single 'Size' column predate the
    # MinSize/MaxSize split and would otherwise raise a bare pandas KeyError.
    # NumMale/NumFemale were renamed to MinNumMale/MinNumFemale, because both were
    # always lower bounds and the old names read as exact counts. Say so, rather than
    # reporting the new names as merely missing.
    renamed = {"NumMale": "MinNumMale", "NumFemale": "MinNumFemale"}
    legacy = [old for old, new in renamed.items()
              if old in corral_defs.columns and new not in corral_defs.columns]
    if legacy:
        raise UsageError(
            "Corral file %s uses %s. Those columns are lower bounds rather than exact "
            "counts, and were renamed to %s to say so. Rename them in the header row; "
            "the values are unchanged. To fix the proportion of each sex rather than "
            "only its floor, see the optional RatioMale/RatioFemale columns."
            % (args.corral_file, " and ".join(legacy),
               " and ".join(renamed[old] for old in legacy)))

    missing_columns = [c for c in CORRAL_COLUMNS if c not in corral_defs.columns]
    if missing_columns:
        raise UsageError(
            "Corral file %s is missing required column(s): %s. Files with a single "
            "'Size' column predate the MinSize/MaxSize split and must be migrated."
            % (args.corral_file, ", ".join(missing_columns)))

    inverted = corral_defs.index[corral_defs['MinSize'] > corral_defs['MaxSize']].tolist()
    if inverted:
        raise UsageError(
            "Corral file %s declares MinSize > MaxSize on row(s) %s, which no "
            "allocation can satisfy." % (args.corral_file, inverted))

    sex_ratios = validate_sex_ratios(corral_defs, args.corral_file)
    species_groups = validate_species_groups(corral_defs, args.corral_file, species,
                                             must_allocate, args.individuals_file)

    LOG.info("Corral definition:\n%s", corral_defs)

    # Size, ordering, non-negativity, symmetry and diagonal, all in one place. Ordering is
    # Ordering is the check that matters most here: it is the one that fails silently.
    connections = validate_relatedness(pr, names, args.pairwise_relatedness_file)

    # objective_function.name, not the member itself: IntEnum.__str__ returns
    # "CobreederObjectiveFunction.ALL_PAIRS" up to Python 3.10 and the bare number "1"
    # from 3.11 on, so printing the member made this field depend on the interpreter.
    # Both supported versions now write the same token, and it is the same token the
    # OBJECTIVE argument accepts. The numeric id remains in the SOLUTION and COMPLETION
    # records as objective_id.
    print("###COBREEDER-ARGS", args.pairwise_relatedness_file, args.corral_file,
          objective_function.name, args.unique_run_id, sep=',')

    return (connections, corral_defs, names, males, females, allocate_first_corral,
            species, allele_counts, must_allocate, objective_function, unique_id,
            sex_ratios, species_groups)


def make_solver(args):
    """A solver configured the same way for bound derivation and the final solve.

    One place to change, so the multi-objective mode's component bounds are never
    derived under different settings from the run that consumes them.

    Reproducibility rests on two settings together, not one. A fixed random_seed alone
    is not enough: with several workers CP-SAT's portfolio communicates as the search
    runs, so timing between threads changes the result. num_workers=1 with a fixed seed
    is bitwise reproducible; num_workers=0 (one per core) is not, whatever the seed.

    A wall-clock limit reintroduces machine dependence, since a faster machine explores
    further before stopping. max_deterministic_time bounds the search by work done
    rather than time elapsed, and is the only limit that keeps a run reproducible.
    """
    solver = cp_model.CpSolver()
    solver.parameters.random_seed = args.seed
    solver.parameters.num_workers = args.workers
    if args.time_limit:
        solver.parameters.max_time_in_seconds = args.time_limit
    if args.deterministic_time:
        solver.parameters.max_deterministic_time = args.deterministic_time
    return solver


def build_model(data, args):
    """Build the CP-SAT model: variables, constraints, and the chosen objective.

    Separated from solving so a model can be constructed and inspected without a
    search, which is what makes the constraints testable in-process. Returns a
    namespace carrying the model and the handles the solve and reporting need.
    """
    (connections, corral_defs, names, males, females, allocate_first_corral,
     species, allele_counts, must_allocate, objective_function, unique_id,
     sex_ratios, species_groups) = data

    num_individuals = len(connections)
    num_corrals = len(corral_defs)
    all_corrals = range(num_corrals)
    all_individuals = range(num_individuals)

    # Create the cp model.
    model = cp_model.CpModel()

    #
    # Decision variables
    #
    placement = {}
    # 1 where corral t will accept individual g at all, and 1 where accepting g would
    # count against t's MaxNumNonComp cap. Both were previously computed and then never
    # referenced by any constraint, which is what made the species columns decorative.
    # A third structure, compulsory_match_individual_corral_violated, is gone: it was
    # unused and tested species[g] != "R", hardcoding one scenario's species code.
    individual_corral_compatibility = {}
    optional_group_allocation = {}
    individual_must_be_allocated = {}
    individual_allele_count = {}
    LOG.debug("Corral, Individual, Species, compatible, counts against MaxNumNonComp")

    for g in all_individuals:
        individual_must_be_allocated[g] = int(must_allocate[g])
        individual_allele_count[g] = allele_counts[g]

    for t in all_corrals:
        preferred, optional = species_groups[t]
        # A corral naming no species accepts every species; see validate_species_groups.
        unrestricted = not preferred and not optional
        for g in all_individuals:
            placement[(t, g)] = model.NewBoolVar("individual %i placed in corral %i" % (g, t))
            compatible = unrestricted or species[g] in preferred or species[g] in optional
            non_compulsory = (not unrestricted) and species[g] in optional
            individual_corral_compatibility[(t, g)] = 1 if compatible else 0
            optional_group_allocation[(t, g)] = 1 if non_compulsory else 0
            LOG.debug("%s %s %s %s %s", t, g, species[g], compatible, non_compulsory)



    colocated = {}
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            colocated[(g1, g2)] = model.NewBoolVar(
                "individual %i colocated with individual %i" % (g1, g2)
            )

    same_corral = {}
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            for t in all_corrals:
                same_corral[(g1, g2, t)] = model.NewBoolVar(
                    "individuals %i and %i are both placed in corral %i" % (g1, g2, t)
                )

    opposing_sex = {}
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            opposing_sex[(g1, g2)] = 0 if males[g1] == males[g2] else 1

    # The two components the multi-objective mode combines are kept as explicit
    # (coefficient, variable) lists as well as expressions, because normalising them
    # means rescaling their coefficients -- see normalise_terms.
    allele_terms = [
        (individual_allele_count[g], placement[(t, g)])
        for g in range(num_individuals)
        for t in range(num_corrals)
    ]
    prsq_terms = [
        (connections[g1][g2] * connections[g1][g2], colocated[g1, g2])
        for g1 in range(num_individuals - 1)
        for g2 in range(g1 + 1, num_individuals)
        if connections[g1][g2] > 0
    ]

    # Everything an objective builder needs. `pairs` is pre-filtered to the pairs with
    # non-zero relatedness, which is a sparsity optimisation only: every objective weights
    # each pair by its relatedness, so a zero coefficient contributes nothing to any of
    # them. Negative values are rejected on load, so nothing meaningful is excluded here.
    # An objective that counted co-located pairs rather than weighting them would need the
    # unfiltered list.
    objective_ctx = SimpleNamespace(
        allele_terms=allele_terms,
        prsq_terms=prsq_terms,
        connections=connections,
        colocated=colocated,
        opposing_sex=opposing_sex,
        pairs=[(g1, g2)
               for g1 in range(num_individuals - 1)
               for g2 in range(g1 + 1, num_individuals)
               if connections[g1][g2] > 0],
    )

    #
    # Constraints
    #

    total_allocated = sum(
        placement[(t, g)]
        for g in range(num_individuals)
        for t in range(num_corrals)
    )
    model.Add(total_allocated >= args.total_individuals)

    # Setting individual-centric constraints.
    for g in all_individuals:

        # All individuals which must be placed are allocated to a corral.
        # Else, each individual is placed in most one corral.
        if individual_must_be_allocated[g]:
            LOG.debug("individual %s must be allocated", g)
            model.Add(sum(placement[(t, g)] for t in all_corrals) == 1)
        else:
            model.Add(sum(placement[(t, g)] for t in all_corrals) <= 1)

    # Setting corral-centric constraints
    for t in all_corrals:
        preferred, optional = species_groups[t]
        LOG.info("Corral %s is of size %i to %i; preferred species %s, optional %s",
                 t, corral_defs['MinSize'][t], corral_defs['MaxSize'][t],
                 ", ".join(sorted(preferred)) or "any",
                 ", ".join(sorted(optional)) or "none")

        # Each corral's occupancy lies within its declared [MinSize, MaxSize] range.
        model.Add(sum(placement[(t, g)] for g in all_individuals) >= corral_defs['MinSize'][t])
        model.Add(sum(placement[(t, g)] for g in all_individuals) <= corral_defs['MaxSize'][t])

        # Each corral holds at least this many males. A floor, not an exact count:
        # use the optional RatioMale/RatioFemale columns to fix the proportion.
        model.Add(sum(males[g] * placement[(t, g)] for g in all_individuals)
                  >= corral_defs['MinNumMale'][t])

        # Each corral holds at least this many females.
        model.Add(sum(females[g] * placement[(t, g)] for g in all_individuals)
                  >= corral_defs['MinNumFemale'][t])

        # Optional exact male:female ratio, for corrals whose RatioMale/RatioFemale
        # are set. Stated as a cross-multiplication so it stays linear: a ratio of
        # m:f means males/females == m/f, hence f*males == m*females.
        #
        # This sits on top of the MinNumMale/MinNumFemale floors rather than replacing
        # them: the minimums set a floor on each sex, the ratio fixes the proportion
        # between them. Combinations that cannot hold together are rejected at load
        # time by validate_sex_ratios, so reaching here means some occupancy satisfies
        # both.
        if t in sex_ratios:
            ratio_male, ratio_female, ratio_mode = sex_ratios[t]
            LOG.info("Corral %s requires a %d:%d male:female ratio (%s); satisfiable "
                     "at occupancy %s", t, ratio_male, ratio_female, ratio_mode,
                     sex_ratio_occupancies(
                                        ratio_male, ratio_female, ratio_mode,
                         int(corral_defs['MinSize'][t]),
                         int(corral_defs['MaxSize'][t]),
                         int(corral_defs['MinNumMale'][t]),
                         int(corral_defs['MinNumFemale'][t])))
            male_side = ratio_female * sum(
                males[g] * placement[(t, g)] for g in all_individuals)
            female_side = ratio_male * sum(
                females[g] * placement[(t, g)] for g in all_individuals)
            if ratio_mode == "exact":
                model.Add(male_side == female_side)
            elif ratio_mode == "min-females":
                model.Add(male_side <= female_side)
            else:                                            # min-males
                model.Add(male_side >= female_side)


        # Species compatibility. An individual whose species is neither preferred nor
        # optional here simply cannot be placed here. Stated as a fixed zero rather than
        # a sum over compatible placements, so presolve removes the variable outright
        # instead of carrying it through the search.
        #
        # This is the constraint that was missing: individual_corral_compatibility was
        # computed and never used, so Species had no effect on any allocation. Every
        # corral configuration shipped in this repository lists all of its cohort's
        # species, so enforcing it changes none of them -- verified in the tests.
        incompatible = [g for g in all_individuals
                        if not individual_corral_compatibility[(t, g)]]
        if incompatible:
            LOG.info("Corral %s excludes %d individual(s) on species grounds", t,
                     len(incompatible))
            for g in incompatible:
                model.Add(placement[(t, g)] == 0)

        # Cap how many occupants come from the optional species rather than the
        # preferred ones. -1 means no cap. Previously commented out.
        cap = corral_defs['MaxNumNonComp'][t]
        if cap != -1:
            optional_here = sum(optional_group_allocation[(t, g)] * placement[(t, g)]
                                for g in all_individuals)
            LOG.info("Corral %s admits at most %d individual(s) of its optional "
                     "species", t, int(cap))
            model.Add(optional_here <= int(cap))

        # Add 'maximum pairwise relatedness' constraint for corrals whose MaxPR value != -1.
        if corral_defs['MaxPR'][t] != -1:
            model.Add(
                sum(
                    connections[g1][g2] * same_corral[(g1, g2, t)]
                    for g1 in range(num_individuals - 1)
                    for g2 in range(g1 + 1, num_individuals)
                    if connections[g1][g2] > corral_defs['MaxPR'][t]
                ) < 1
            )

    # Link colocated with placement decisions
    for g1 in range(num_individuals - 1):
        for g2 in range(g1 + 1, num_individuals):
            for t in all_corrals:
                # Link same_corral and placement.
                model.AddBoolOr(
                    [
                        placement[(t, g1)].Not(),
                        placement[(t, g2)].Not(),
                        same_corral[(g1, g2, t)],
                    ]
                )
                model.AddImplication(same_corral[(g1, g2, t)], placement[(t, g1)])
                model.AddImplication(same_corral[(g1, g2, t)], placement[(t, g2)])

            # Link colocated and same_corral.
            model.Add(
                sum(same_corral[(g1, g2, t)] for t in all_corrals) == colocated[(g1, g2)]
            )

    # Individuals carrying a pre-assigned corral are pinned to it. This is the only
    # place that happens; a second, identical copy inside the individual loop above
    # was removed. Note it is not symmetry breaking, which the original comment here
    # claimed -- corral symmetry is not broken anywhere.
    pinned = [g for g in range(len(allocate_first_corral))
              if allocate_first_corral[g] != -1]
    if pinned:
        # Named rather than counted, and the semantics restated, because the column is
        # called AssignToFirstCorral while its value is a zero-based corral index: 0 pins
        # to the first corral and 1 pins to the second. Four shipped cohorts set it, two
        # of them to 1 -- see the pre-allocation note in the README.
        LOG.info("Applying %d pre-assigned corral allocation(s). AssignToFirstCorral is "
                 "a zero-based corral index, not a flag.", len(pinned))
        for g1 in pinned:
            LOG.info("  individual %i (%s) pinned to corral index %i", g1, names[g1],
                     allocate_first_corral[g1])
            model.Add(placement[(allocate_first_corral[g1], g1)] == 1)


    #
    # Objective
    #
    # Set after every constraint is in place. The multi-objective mode derives each
    # component's achievable range by optimising it alone on this model, so the model
    # must already be fully constrained or those ranges would not be achievable.
    #

    if objective_function in OBJECTIVE_BUILDERS:
        sense, builder = OBJECTIVE_BUILDERS[objective_function]
        expression = builder(objective_ctx)
        (model.Maximize if sense == "max" else model.Minimize)(expression)
    elif objective_function == CobreederObjectiveFunction.WEIGHTED_ALLELES_PR_50_50:
        component_terms = {"alleles": allele_terms,
                           "all_pairs_pr_squared": prsq_terms}
        component_exprs = {name: build(objective_ctx)
                           for name, build in WEIGHTED_COMPONENT_BUILDERS.items()}
        weights = {"alleles": args.weight_a,
                   "all_pairs_pr_squared": args.weight_b}

        if args.objective_bounds == "auto":
            LOG.info("Deriving component bounds from the data (%i solves; use "
                     "--objective-bounds to supply them instead).",
                     2 * len(component_exprs))
            bounds = derive_component_bounds(model, lambda: make_solver(args),
                                             component_exprs,
                                             args.bounds_time_limit)
        else:
            try:
                values = [int(v) for v in args.objective_bounds.split(",")]
            except ValueError:
                raise UsageError(
                    "--objective-bounds must be 'auto' or four integers "
                    "LO_A,HI_A,LO_B,HI_B; got %r" % args.objective_bounds)
            if len(values) != 4:
                raise UsageError(
                    "--objective-bounds needs exactly four integers "
                    "LO_A,HI_A,LO_B,HI_B; got %i" % len(values))
            bounds = {("alleles", "min"): values[0], ("alleles", "max"): values[1],
                      ("all_pairs_pr_squared", "min"): values[2],
                      ("all_pairs_pr_squared", "max"): values[3]}

        objective = 0
        total_magnitude = 0
        for name, direction in WEIGHTED_COMPONENTS:
            lo, hi = bounds[(name, "min")], bounds[(name, "max")]
            scaled, offset, magnitude = normalise_terms(
                component_terms[name], lo, hi, weights[name],
                args.objective_scale, name)
            objective += direction * (sum(c * v for c, v in scaled) - offset)
            total_magnitude += magnitude
            print("###COBREEDER-BOUNDS,%s,%d,%d,%d,%d,%s"
                  % (name, lo, hi, weights[name], args.objective_scale,
                     args.unique_run_id))

        # Fail with an explanation rather than letting CP-SAT reject the model.
        if total_magnitude > INT64_MAX:
            raise UsageError(
                "The weighted objective could reach %.3e, above the int64 limit of "
                "%.3e. Lower --objective-scale (currently %d) or the weights."
                % (total_magnitude, INT64_MAX, args.objective_scale))

        LOG.info("Weighted objective magnitude bound %.3e (%.0fx headroom under "
                 "int64).", total_magnitude, INT64_MAX / max(total_magnitude, 1))
        model.Maximize(objective)


    return SimpleNamespace(
        model=model,
        placement=placement,
        colocated=colocated,
        same_corral=same_corral,
        names=names,
        males=males,
        females=females,
        sex_ratios=sex_ratios,
        num_individuals=num_individuals,
        num_corrals=num_corrals,
        objective_function=objective_function,
        unique_id=unique_id,
    )


def run_solver(bundle, args):
    """Search the model, streaming intermediate solutions through the printer."""
    num_corrals = bundle.num_corrals
    num_individuals = bundle.num_individuals
    objective_function = bundle.objective_function
    model = bundle.model
    placement = bundle.placement
    names = bundle.names
    unique_id = bundle.unique_id

    print("###COBREEDER-SEARCH,%i,%i,%f,%f,%s"
          % (args.seed, args.workers, args.time_limit, args.deterministic_time,
             args.unique_run_id))

    paramstring = "%i,%i,%i" % (num_corrals, objective_function, num_individuals)

    # Solve model.
    solver = make_solver(args)
    solution_printer = CobreederPrinter(placement, names, num_corrals, num_individuals, paramstring, unique_id)

    # CP-SAT's search log is the single largest contributor to output volume, and it
    # went to stdout alongside the ### records. Route it through the logger instead,
    # and only ask for it when someone is listening.
    if LOG.isEnabledFor(logging.DEBUG):
        solver.parameters.log_search_progress = True
        solver.parameters.log_to_stdout = False
        solver.log_callback = lambda line: LOG.debug("cp-sat: %s", line.rstrip())

    status = solver.Solve(model, solution_printer)


    return status, solver, solution_printer, paramstring


def report_outcome(status, solver, solution_printer, paramstring, bundle, args):
    """Print the statistics and completion record; return the process exit code."""
    males = bundle.males
    females = bundle.females
    sex_ratios = bundle.sex_ratios
    num_corrals = bundle.num_corrals

    status_name, exit_code = SOLVER_STATUS.get(
        status, ("UNRECOGNISED_STATUS_%i" % status, 5))

    print("Statistics")
    print("  - status       : %s (%i)" % (status_name, status))
    print("  - conflicts    : %i" % solver.NumConflicts())
    print("  - branches     : %i" % solver.NumBranches())
    print("  - wall time    : %f s" % solver.WallTime())
    print("  - num solutions: %i" % solution_printer.num_solutions())

    # The completion record keeps the numeric status as its first field; consumers
    # parse it positionally, so the layout must not change.
    print("###COBREEDER-COMPLETION,%i,%i,%i,%f,%i,%s,%s" % (status, solver.NumConflicts(),
                                                        solver.NumBranches(),
                                                        solver.WallTime(),
                                                        solution_printer.num_solutions(),
                                                        paramstring,
                                                        args.unique_run_id
                                                        ))

    # Say which of the three no-solution cases occurred. Previously all of them
    # printed "No solution found.", so a rejected model, a proven-infeasible
    # instance and an inconclusive search were indistinguishable.
    if status in (cp_model.OPTIMAL, cp_model.FEASIBLE):
        print("%s: %i improving solution(s) found." % (status_name,
                                                      solution_printer.num_solutions()))
    elif status == cp_model.INFEASIBLE:
        print("INFEASIBLE: no allocation satisfies these constraints. This is a "
              "definitive result about the inputs, not a failure.")
        if sex_ratios:
            cohort_males, cohort_females = sum(males), sum(females)
            print("  Sex ratios are set on %d of %d corrals. The cohort holds %d "
                  "male and %d female individuals."
                  % (len(sex_ratios), num_corrals, cohort_males, cohort_females))
            cap = sex_ratio_placement_cap(sex_ratios, num_corrals,
                                          cohort_males, cohort_females)
            if cap is not None:
                print("  With this ratio on every corral, at most %d individuals can "
                      "be placed whatever the corral sizes, because the scarcer sex "
                      "runs out first. total_individuals is %d."
                      % (cap, args.total_individuals))
                if args.total_individuals > cap:
                    print("  That alone makes the model infeasible: lower "
                          "total_individuals to %d or less, or change the ratio."
                          % cap)
    elif status == cp_model.MODEL_INVALID:
        print("MODEL_INVALID: the model was rejected before solving, so this run "
              "produced no result. See the solver log above; an integer overflow "
              "in the objective is the usual cause.", file=sys.stderr)
    elif status == cp_model.UNKNOWN:
        print("UNKNOWN: the search ended without proving optimality or "
              "infeasibility.", file=sys.stderr)
    else:
        print("%s: unrecognised solver status." % status_name, file=sys.stderr)

    return exit_code


def solve_with_discrete_model(args):
    """Load, build, solve, report. Returns the exit code; see SOLVER_STATUS."""
    data = build_data(args)
    bundle = build_model(data, args)
    status, solver, solution_printer, paramstring = run_solver(bundle, args)
    return report_outcome(status, solver, solution_printer, paramstring, bundle, args)

def objective_help():
    """One line per objective function, for the epilog.

    Built from the enum rather than typed out, so a new objective cannot be added
    without appearing here.
    """
    summaries = {
        # Both of these minimise summed relatedness; the summaries used to read
        # "maximise co-located pairs", which named the wrong quantity and the wrong
        # direction. OBJECTIVE_BUILDERS registers each as ("min", ...), and
        # test_help_summaries_name_the_direction_the_registry_uses now checks the
        # summaries against the registry so the two cannot drift apart again.
        "ALL_PAIRS": "minimise summed relatedness over co-located pairs",
        "MALE_FEMALE": "minimise summed relatedness over opposite-sex co-located "
                       "pairs",
        "ALLELES_MIN": "minimise total ghost alleles housed",
        "ALLELES_MAX": "maximise total ghost alleles housed",
        "ALL_PAIRS_PR_MIN": "minimise summed relatedness over co-located pairs",
        "ALL_PAIRS_PR_MAX": "maximise summed relatedness over co-located pairs",
        "ALL_PAIRS_PR_MIN_SQUARED": "as ALL_PAIRS_PR_MIN, squared -- penalises close kin "
                                    "disproportionately",
        "ALL_PAIRS_PR_MAX_SQUARED": "as ALL_PAIRS_PR_MAX, squared",
        "MALE_FEMALE_SQUARED": "minimise squared relatedness over opposite-sex pairs",
        "WEIGHTED_ALLELES_PR_50_50": "weighted sum of alleles and relatedness; needs "
                                     "weight_a and weight_b",
    }
    width = max(len(e.name) for e in CobreederObjectiveFunction)
    return "\n".join(
        "  %-*s  %s" % (width, e.name, summaries.get(e.name, ""))
        for e in CobreederObjectiveFunction)


def epilog(subcommand=""):
    """Text for the --help footer.

    subcommand is "run " for the top-level parser and "" for the run subparser, because
    argparse sets %(prog)s to "solver.py" in the first case and "solver.py run" in the
    second -- so the examples need the command name in one and not the other.

    argparse expands the epilog with `text % dict(prog=...)`, so every literal percent
    sign here has to be doubled. Building this by concatenation rather than %-formatting
    keeps that escaping in one place: formatting the objective list in with % would consume
    one level of escaping and make --help raise TypeError.
    """
    return (
        "objective functions:\n"
        + objective_help() + "\n"
        "\n"
        "examples:\n"
        "  Single objective, the canid case study:\n"
        "    %(prog)s " + subcommand + "data/wolves/individuals.csv data/wolves/pr-scaled.csv \\\n"
        "        data/wolves/corrals-a.csv ALL_PAIRS_PR_MIN_SQUARED demo 0 0 10\n"
        "\n"
        "  Weighted objective, 70%% alleles against 30%% relatedness:\n"
        "    %(prog)s " + subcommand + "data/wolves/individuals.csv data/wolves/pr-scaled.csv \\\n"
        "        data/wolves/corrals-a.csv WEIGHTED_ALLELES_PR_50_50 demo 70 30 10\n"
        "\n"
        "  Reproducible across machines, bounded by work rather than wall-clock:\n"
        "    %(prog)s " + subcommand + "... --seed 1 --workers 1 --deterministic-time 60\n"
        "\n"
        "Results go to stdout as ###COBREEDER-* records; developer traces go to stderr.\n"
        "Exit status carries the solver status: 0 for any definitive answer including\n"
        "INFEASIBLE, 1 bad input, 2 command-line misuse, 3 MODEL_INVALID, 4 UNKNOWN.\n")


def build_arg_parser():
    """The command-line interface.

    prog is left to argparse, which derives it from argv[0] -- correct both for
    `python3 cobreeder/solver.py` and for the installed `cobreeder` script. It was
    previously the literal placeholder 'PROG'.
    """
    parser = argparse.ArgumentParser(
        description="Design group-living captive breeding programmes by constraint\n"
                    "optimisation. Allocates individuals to breeding groups subject to\n"
                    "group size, sex composition and relatedness constraints.",
        epilog=epilog(subcommand="run "),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--version", action="version",
                       version="cobreeder %s" % __version__)

    # required=True so that no arguments produces a usage message and exit 2. Without
    # it, argparse returned an empty namespace and the crash surfaced from inside
    # build_data as an AttributeError -- at exit 0, reporting success for a run that
    # never happened.
    subparsers = parser.add_subparsers(dest="command", metavar="COMMAND", required=True,
                                       title="commands")
    run_parser = subparsers.add_parser(
        "run", help="solve one allocation problem",
        description="Solve one allocation problem and report every improving solution.",
        epilog=epilog(), formatter_class=argparse.RawDescriptionHelpFormatter)

    inputs = run_parser.add_argument_group(
        "inputs",
        "All three are CSV. See the README for each schema; the relatedness matrix\n"
        "is validated against the individuals file on load.")
    inputs.add_argument("individuals_file", metavar="INDIVIDUALS",
                        help="the cohort: Name, Male, Female, AssignToFirstCorral, "
                             "Species, Alleles, MustAllocate")
    inputs.add_argument("pairwise_relatedness_file", metavar="RELATEDNESS",
                        help="square, symmetric, non-negative integer relatedness "
                             "matrix, in the same order as INDIVIDUALS")
    inputs.add_argument("corral_file", metavar="CORRALS",
                        help="one row per breeding group: sizes, sex floors, MaxPR")

    problem = run_parser.add_argument_group("problem")
    problem.add_argument("obj_function", metavar="OBJECTIVE",
                         choices=[e.name for e in CobreederObjectiveFunction],
                         help="what to optimise; one of the names listed below")
    problem.add_argument("unique_run_id", metavar="RUN_ID",
                         help="label stamped into every output record, so results from "
                              "several runs can be concatenated and told apart")
    problem.add_argument("weight_a", type=int, metavar="WEIGHT_A",
                         help="weight on the ghost-allele component. Used only by "
                              "WEIGHTED_ALLELES_PR_50_50; pass 0 otherwise. Only the "
                              "ratio to WEIGHT_B matters, so 70 30 and 7 3 are the same "
                              "run")
    problem.add_argument("weight_b", type=int, metavar="WEIGHT_B",
                         help="weight on the relatedness component. Used only by "
                              "WEIGHTED_ALLELES_PR_50_50; pass 0 otherwise")
    problem.add_argument("total_individuals", type=int, metavar="MIN_TOTAL",
                         help="minimum individuals allocated across all groups. Group "
                              "size constraints usually imply this already; it binds "
                              "when groups can hold a range")

    search = run_parser.add_argument_group(
        "search and reproducibility",
        "A fixed seed alone is not enough: with several workers CP-SAT's portfolio\n"
        "communicates as the search runs, so thread timing changes the answer.")
    search.add_argument("--seed", type=int, default=0, metavar="N",
                        help="CP-SAT random seed. Default 0.")
    search.add_argument("--workers", type=int, default=1, metavar="N",
                        help="Search workers. Default 1, which with a fixed seed makes "
                             "runs bitwise reproducible. 0 uses one worker per core, "
                             "which is faster on hard instances but gives a different "
                             "answer each time.")
    search.add_argument("--time-limit", type=float, default=0.0, metavar="SECONDS",
                        help="Wall-clock limit; 0 means none. Note this reintroduces "
                             "machine dependence, because a faster machine searches "
                             "further before stopping.")
    search.add_argument("--deterministic-time", type=float, default=0.0, metavar="UNITS",
                        help="Limit the search by work done rather than time elapsed; 0 "
                             "means none. The only limit that keeps a run reproducible "
                             "across machines.")

    output = run_parser.add_argument_group("output")
    output.add_argument("--verbose", "-v", action="store_true",
                        help="Add developer traces to stderr: the per-corral "
                             "compatibility dump, the raw linear expressions, and "
                             "CP-SAT's search log. All are quadratic in population "
                             "size, so they are off by default.")
    output.add_argument("--quiet", "-q", action="store_true",
                        help="Suppress the stderr traces entirely, leaving only "
                             "warnings. stdout is unaffected either way.")

    weighted = run_parser.add_argument_group(
        "weighted objective", "Affect WEIGHTED_ALLELES_PR_50_50 only.")
    weighted.add_argument("--objective-bounds", type=str, default="auto",
                          metavar="auto|LO_A,HI_A,LO_B,HI_B",
                          help="Component ranges used to normalise the weighted "
                               "objective. 'auto' (default) derives them from the data "
                               "by optimising each component alone, so no "
                               "scenario-specific value is baked into the code. Supply "
                               "four integers to skip those solves.")
    weighted.add_argument("--objective-scale", type=int, default=1000000, metavar="N",
                          help="Integer span each normalised component is scaled to. "
                               "Larger is more precise, smaller leaves more headroom "
                               "under the int64 objective limit. Default 1000000.")
    weighted.add_argument("--bounds-time-limit", type=float, default=0.0,
                          metavar="SECONDS",
                          help="Per-solve time limit when deriving bounds. Each must "
                               "still prove optimality, so this only bounds how long a "
                               "failure takes. 0 means no limit.")
    return parser


def main(argv: Sequence[str]) -> int:
    parser = build_arg_parser()
    args = parser.parse_args(argv[1:])
    configure_logging(verbose=getattr(args, "verbose", False),
                      quiet=getattr(args, "quiet", False))

    # Returned to __main__ and passed to sys.exit, so the solver status reaches
    # the shell. See SOLVER_STATUS for the mapping.
    return solve_with_discrete_model(args)


def cli():
    """Console-script entry point, installed as `cobreeder` by pyproject.toml.

    Takes no arguments so Poetry can wire it directly, and exits with the solver
    status from SOLVER_STATUS. Input problems are reported as a message rather than
    a traceback.
    """
    # There was a sys.setrecursionlimit(0x100000) here, paired with a commented-out
    # resource.setrlimit call annotated "Will segfault without this line." Neither was
    # needed. Model building does not recurse in proportion to problem size: the smallest
    # workable limit is 21-22 for every cohort from 10 to 51 individuals, 45 to 1275
    # pairs, and a complete run -- build, solve, per-solution callback, reporting -- needs
    # 37, again independent of size. Python's default of 1000 is 27 times that, so the
    # line raised the ceiling roughly 28,000-fold to no purpose.
    # tests/test_weighted_objective.py::TestRecursionDepth keeps that true.
    try:
        sys.exit(main(sys.argv))
    except UsageError as err:
        print("error: %s" % err, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    cli()
