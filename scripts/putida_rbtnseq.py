"""
Helper API for the P. putida KT2440 RB-TnSeq skill.

Loads from a parquet cache built by ``build_cache.py``. Each loader is memoized —
first call is ~50 ms, subsequent calls are free.

Data location
-------------
Resolved in this order:

1. ``$PUTIDA_RBTNSEQ_DATA``  (environment variable)
2. ``<repo root>/data``      (default; see docs/DATA.md for how to populate it)

Matrices
--------
Canonical fitness matrix = the **paper** version (332 experiments, Borchert 2024).
The public Fitness Browser snapshot (314 experiments) is available on the matrix
functions via ``use_paper=False``.

Cofitness comes from the Fitness Browser's precomputed table (314-experiment
basis); the paper does not ship a recomputed cofitness table.

Standard significance call: ``|fit| > 1 AND |t| > 4``.

Public API
----------
    gene(locus)                         -> Series of annotation + reanno + UniProt
    fitness(locus)                      -> Series indexed by expName
    t_stat(locus)                       -> Series indexed by expName
    phenotypes(locus, ...)              -> DataFrame of significant conditions
    cofitness(locus, n)                 -> DataFrame of top cofit neighbors
    experiments(group, contains)        -> DataFrame of filtered metadata
    condition_top_genes(expName, ...)   -> DataFrame of top phenotypes in a condition
    condition_compare(a, b)             -> DataFrame of fit_a, fit_b, delta
    substrate_grid(genes, substrates)   -> DataFrame (genes x experiments)
    qc_pass(expName)                    -> (bool, dict of metrics)
    has_signed_t(expName)               -> bool (see SETS_WITHOUT_SIGNED_T)
    plot_gene_profile(locus, out=...)   -> PNG of a gene's fitness across experiments
    plot_substrate_heatmap(grid, out=)  -> PNG heatmap of a substrate grid
"""

from __future__ import annotations

import os
from functools import lru_cache
from pathlib import Path
from typing import Iterable

import pandas as pd

# --------------------------------------------------------------------------
# data location
# --------------------------------------------------------------------------

_REPO_ROOT = Path(__file__).resolve().parents[1]


def data_dir() -> Path:
    """Directory holding the raw source files and the ``cache/`` subdirectory."""
    env = os.environ.get("PUTIDA_RBTNSEQ_DATA")
    return Path(env).expanduser() if env else _REPO_ROOT / "data"


def cache_dir() -> Path:
    return data_dir() / "cache"


class DataNotFound(FileNotFoundError):
    """Raised when the parquet cache is missing, with instructions to build it."""


def _read(name: str) -> pd.DataFrame:
    path = cache_dir() / name
    if not path.exists():
        raise DataNotFound(
            f"Missing cache file: {path}\n\n"
            "The dataset is not redistributed with this repository.\n"
            "  1. Download the source files  -> see docs/DATA.md\n"
            "  2. Build the cache            -> python scripts/build_cache.py\n"
            "  3. Or point at an existing copy:\n"
            "       export PUTIDA_RBTNSEQ_DATA=/path/to/rbtnseq"
        )
    return pd.read_parquet(path)


# Significance defaults (Arkin lab convention)
FIT_THRESH = 1.0
T_THRESH = 4.0

# Experiment sets whose `t` values are stored as UNSIGNED magnitudes in the paper
# xlsx (verified: 0 negative values across 369,096 measurements; every other set
# is ~52% negative). Use |t| for confidence and signed `fit` for direction.
# See docs/DATA_CAVEATS.md.
SETS_WITHOUT_SIGNED_T = frozenset({"set100", "set101"})


# ---------------------------------------------------------------------------
# loaders (memoized — these run once per process)
# ---------------------------------------------------------------------------

@lru_cache(maxsize=1)
def _genes() -> pd.DataFrame:
    return _read("genes.parquet")


@lru_cache(maxsize=1)
def _reanno() -> pd.DataFrame:
    return _read("reanno.parquet")


@lru_cache(maxsize=1)
def _uniprot() -> pd.DataFrame:
    return _read("uniprot.parquet")


@lru_cache(maxsize=1)
def _metadata_paper() -> pd.DataFrame:
    return _read("metadata_paper.parquet")


@lru_cache(maxsize=1)
def _metadata_public() -> pd.DataFrame:
    return _read("metadata.parquet")


@lru_cache(maxsize=1)
def _fit_paper() -> pd.DataFrame:
    return _read("fit_paper.parquet")


@lru_cache(maxsize=1)
def _t_paper() -> pd.DataFrame:
    return _read("t_paper.parquet")


@lru_cache(maxsize=1)
def _long_paper() -> pd.DataFrame:
    return _read("fit_paper_long.parquet")


@lru_cache(maxsize=1)
def _fit_public() -> pd.DataFrame:
    return _read("fit.parquet")


@lru_cache(maxsize=1)
def _t_public() -> pd.DataFrame:
    return _read("t.parquet")


@lru_cache(maxsize=1)
def _long_public() -> pd.DataFrame:
    return _read("fit_long.parquet")


@lru_cache(maxsize=1)
def _cofit() -> pd.DataFrame:
    return _read("cofit.parquet")


@lru_cache(maxsize=1)
def _specific() -> pd.DataFrame:
    return _read("specific_phenotypes.parquet")


def _fit(use_paper: bool = True) -> pd.DataFrame:
    return _fit_paper() if use_paper else _fit_public()


def _tmat(use_paper: bool = True) -> pd.DataFrame:
    return _t_paper() if use_paper else _t_public()


def _long(use_paper: bool = True) -> pd.DataFrame:
    return _long_paper() if use_paper else _long_public()


def _metadata(use_paper: bool = True) -> pd.DataFrame:
    return _metadata_paper() if use_paper else _metadata_public()


META_COLS = ["orgId", "locusId", "sysName", "geneName", "desc"]

PHENOTYPE_COLS = ["expName", "expDesc", "expGroup", "fit", "t"]
PHENOTYPE_DTYPES = {
    "expName": "string", "expDesc": "string", "expGroup": "string",
    "fit": "float64", "t": "float64",
}


def _empty_phenotypes() -> pd.DataFrame:
    """Empty result with the same columns *and dtypes* as a populated one, so
    results concatenate cleanly across genes."""
    return pd.DataFrame(
        {c: pd.Series(dtype=PHENOTYPE_DTYPES[c]) for c in PHENOTYPE_COLS}
    )


# ---------------------------------------------------------------------------
# gene-level
# ---------------------------------------------------------------------------

def gene(locus: str) -> pd.Series:
    """Full annotation for one gene: gene table + reanno + UniProt + fitness summary.

    Genes present in the genome annotation but absent from the fitness matrix
    (essential genes, or genes with too few transposon insertions to score —
    929 of 5,661 in this dataset) return annotation with ``has_fitness_data=False``
    and null fitness summary fields, rather than raising.

    >>> gene("PP_0154")["desc"]
    'propionyl-CoA:succinate CoA transferase'
    """
    row = _genes().loc[_genes()["locusId"] == locus]
    if row.empty:
        raise KeyError(f"locusId not found in gene table: {locus}")
    out = row.iloc[0].to_dict()

    re_row = _reanno().loc[_reanno()["locusId"] == locus]
    if not re_row.empty:
        out["reannotation"] = re_row.iloc[0]["reannotation"]
        out["reanno_comment"] = re_row.iloc[0]["comment"]

    up = _uniprot().loc[_uniprot()["locusId"] == locus]
    if not up.empty:
        for col in [
            "Entry", "Protein names", "EC number", "Pathway", "Function [CC]",
            "Gene Ontology (biological process)",
            "Gene Ontology (molecular function)",
            "Gene Ontology (cellular component)",
            "Subcellular location [CC]",
        ]:
            if col in up.columns:
                out[col] = up.iloc[0][col]

    # fitness summary — absent for unscored genes
    try:
        f = fitness(locus)
    except KeyError:
        out["has_fitness_data"] = False
        out["n_experiments"] = 0
        for k in ("mean_fit", "min_fit", "max_fit"):
            out[k] = float("nan")
        out["n_phenotypes"] = 0
    else:
        out["has_fitness_data"] = True
        out["n_experiments"] = int(f.notna().sum())
        out["mean_fit"] = float(f.mean())
        out["min_fit"] = float(f.min())
        out["max_fit"] = float(f.max())
        out["n_phenotypes"] = int(len(phenotypes(locus)))

    return pd.Series(out)


def fitness(locus: str, use_paper: bool = True) -> pd.Series:
    """Fitness vector for one gene, indexed by expName."""
    long = _long(use_paper)
    sub = long.loc[long["locusId"] == locus, ["expName", "fit"]]
    if sub.empty:
        raise KeyError(f"locusId not found in fitness matrix: {locus}")
    return sub.set_index("expName")["fit"].astype(float)


def t_stat(locus: str, use_paper: bool = True) -> pd.Series:
    """t-statistic vector for one gene, indexed by expName.

    Note: values from sets in ``SETS_WITHOUT_SIGNED_T`` are unsigned magnitudes.
    """
    long = _long(use_paper)
    sub = long.loc[long["locusId"] == locus, ["expName", "t"]]
    if sub.empty:
        raise KeyError(f"locusId not found in fitness matrix: {locus}")
    return sub.set_index("expName")["t"].astype(float)


def phenotypes(
    locus: str,
    fit_thresh: float = FIT_THRESH,
    t_thresh: float = T_THRESH,
    use_paper: bool = True,
) -> pd.DataFrame:
    """Significant conditions for a gene (|fit| > fit_thresh AND |t| > t_thresh).

    Returns a DataFrame with columns ``PHENOTYPE_COLS``, sorted by signed fit
    (most negative first). Column order is identical whether or not any rows
    are returned, so results can be concatenated across genes.
    """
    long = _long(use_paper)
    sub = long.loc[long["locusId"] == locus]
    sig = sub.loc[(sub["fit"].abs() > fit_thresh) & (sub["t"].abs() > t_thresh)].copy()
    if sig.empty:
        return _empty_phenotypes()
    meta = _metadata(use_paper)[["expName", "expGroup"]]
    sig = sig.merge(meta, on="expName", how="left")
    return sig[PHENOTYPE_COLS].sort_values("fit").reset_index(drop=True)


# ---------------------------------------------------------------------------
# cofitness
# ---------------------------------------------------------------------------

def cofitness(locus: str, n: int = 20) -> pd.DataFrame:
    """Top-n cofit neighbors of a gene (precomputed by the Fitness Browser).

    Returns a DataFrame with hitId, hitSysName, hitDesc, cofit, rank, conserved,
    sorted by rank ascending (strongest first).
    """
    co = _cofit()
    sub = co.loc[co["locusId"] == locus].sort_values("rank").head(n)
    return sub[["hitId", "hitSysName", "hitDesc", "cofit", "rank", "conserved"]] \
        .reset_index(drop=True)


# ---------------------------------------------------------------------------
# experiments
# ---------------------------------------------------------------------------

def has_signed_t(expName: str) -> bool:
    """False for experiments whose `t` is stored as an unsigned magnitude.

    Affects 78 of the 332 paper experiments (all of set100 and set101). For those,
    use signed ``fit`` for direction and ``|t|`` only for confidence.
    """
    return not any(expName.startswith(s) for s in SETS_WITHOUT_SIGNED_T)


def experiments(
    group: str | None = None,
    contains: str | None = None,
    use_paper: bool = True,
) -> pd.DataFrame:
    """Filter experiment metadata.

    group:    exact match on expGroup (e.g. "carbon source", "nitrogen source", "stress")
    contains: case-insensitive substring match on expDesc OR condition_1
    """
    df = _metadata(use_paper).copy()
    if group is not None:
        df = df.loc[df["expGroup"] == group]
    if contains is not None:
        needle = contains.lower()
        d1 = df["expDesc"].astype(str).str.lower().str.contains(needle, na=False)
        d2 = df.get("condition_1", pd.Series([""] * len(df), index=df.index)) \
            .astype(str).str.lower().str.contains(needle, na=False)
        df = df.loc[d1 | d2]
    return df.reset_index(drop=True)


def condition_top_genes(
    expName: str,
    n: int = 20,
    sign: str = "neg",
    fit_thresh: float = FIT_THRESH,
    t_thresh: float = T_THRESH,
    use_paper: bool = True,
) -> pd.DataFrame:
    """Top genes by fitness in a single condition.

    sign='neg' -> most depleted (gene important for growth, fit < -fit_thresh)
    sign='pos' -> most enriched (gene loss helped, fit > +fit_thresh)
    sign='abs' -> either, ranked by |fit|

    Returns a DataFrame with locusId, sysName, desc, fit, t.
    """
    long = _long(use_paper)
    sub = long.loc[long["expName"] == expName].copy()
    if sub.empty:
        raise KeyError(f"expName not found: {expName}")
    sub = sub.loc[sub["t"].abs() > t_thresh]
    if sign == "neg":
        sub = sub.loc[sub["fit"] < -fit_thresh].sort_values("fit")
    elif sign == "pos":
        sub = sub.loc[sub["fit"] > fit_thresh].sort_values("fit", ascending=False)
    elif sign == "abs":
        sub = sub.loc[sub["fit"].abs() > fit_thresh] \
            .assign(_abs=lambda d: d["fit"].abs()) \
            .sort_values("_abs", ascending=False) \
            .drop(columns="_abs")
    else:
        raise ValueError("sign must be 'neg', 'pos', or 'abs'")
    desc = _genes()[["locusId", "desc"]]
    return sub.merge(desc, on="locusId", how="left") \
        .head(n)[["locusId", "sysName", "desc", "fit", "t"]] \
        .reset_index(drop=True)


def condition_compare(
    exp_a: str,
    exp_b: str,
    fit_thresh: float = FIT_THRESH,
    t_thresh: float = T_THRESH,
    use_paper: bool = True,
) -> pd.DataFrame:
    """Per-gene comparison of two conditions.

    Returns [locusId, sysName, desc, fit_a, t_a, fit_b, t_b, delta] sorted by
    |delta| descending, where ``delta = fit_a - fit_b``. A negative delta means
    the mutant was more impaired in exp_a (that gene matters more in condition a).
    Pre-filtered to genes significant in at least one of the two conditions.
    """
    long = _long(use_paper)
    a = long.loc[long["expName"] == exp_a, ["locusId", "fit", "t"]].rename(
        columns={"fit": "fit_a", "t": "t_a"}
    )
    b = long.loc[long["expName"] == exp_b, ["locusId", "fit", "t"]].rename(
        columns={"fit": "fit_b", "t": "t_b"}
    )
    if a.empty:
        raise KeyError(f"expName not found: {exp_a}")
    if b.empty:
        raise KeyError(f"expName not found: {exp_b}")
    j = a.merge(b, on="locusId", how="inner")
    sig = (
        ((j["fit_a"].abs() > fit_thresh) & (j["t_a"].abs() > t_thresh))
        | ((j["fit_b"].abs() > fit_thresh) & (j["t_b"].abs() > t_thresh))
    )
    j = j.loc[sig].copy()
    j["delta"] = j["fit_a"] - j["fit_b"]
    desc = _genes()[["locusId", "desc"]]
    sysn = long[["locusId", "sysName"]].drop_duplicates(subset="locusId")
    j = j.merge(desc, on="locusId", how="left").merge(sysn, on="locusId", how="left")
    j = j.loc[j["delta"].abs().sort_values(ascending=False).index]
    return j[["locusId", "sysName", "desc", "fit_a", "t_a", "fit_b", "t_b", "delta"]] \
        .reset_index(drop=True)


# ---------------------------------------------------------------------------
# matrices
# ---------------------------------------------------------------------------

def substrate_grid(
    genes: Iterable[str] | None = None,
    substrates: Iterable[str] | None = None,
    use_paper: bool = True,
) -> pd.DataFrame:
    """Slice the fitness matrix to a [genes x experiments] DataFrame.

    genes:      list of locusIds (default: all)
    substrates: list of case-insensitive substrings matched against the column
                headers (default: all)

    The result is indexed by locusId. Experiment columns keep the source header
    format ``"<expName> <expDesc>"`` (e.g. ``"set1IT084 Acetate Carbon Source (20mM)"``);
    the remaining ``META_COLS`` (orgId, sysName, geneName, desc) are carried through
    as leading columns.
    """
    fit = _fit(use_paper)
    cols = [c for c in fit.columns if c not in META_COLS]
    if substrates is not None:
        needles = [s.lower() for s in substrates]
        cols = [c for c in cols if any(n in c.lower() for n in needles)]
    if not cols:
        raise ValueError("No experiments matched the substrate filter")

    df = fit[META_COLS + cols].copy()
    if genes is not None:
        df = df.loc[df["locusId"].isin(list(genes))]
    df = df.set_index("locusId")
    df.index.name = "locusId"
    return df


# ---------------------------------------------------------------------------
# QC
# ---------------------------------------------------------------------------

# Published criteria for a usable RB-TnSeq experiment (Wetmore et al. 2015, mBio;
# also stated in the Fitness Browser help). Each entry: (threshold, comparison, meaning).
#   "ge"  -> value >= threshold
#   "le"  -> value <= threshold
#   "abs_le" -> |value| <= threshold
QC_METRICS: dict[str, tuple[float, str, str]] = {
    "gMed":   (50.0, "ge",     "median reads per gene (sequencing depth)"),
    "mad12":  (0.5,  "le",     "median absolute difference in fitness between gene halves"),
    "cor12":  (0.1,  "ge",     "correlation in fitness between the two halves of each gene"),
    "gccor":  (0.2,  "abs_le", "correlation between gene GC content and fitness (bias check)"),
    "adjcor": (0.25, "abs_le", "correlation between adjacent genes (bias check)"),
}


def qc_pass(expName: str) -> tuple[bool, dict]:
    """Did this experiment pass the published QC thresholds?

    Returns ``(overall_pass, {metric: (value, threshold, passed)})``.

    Uses the **public snapshot** metadata — the paper xlsx does not ship QC
    fields, so the 78 paper-only experiments (set100/set101) cannot be QC'd and
    raise KeyError.

    Note: the public Fitness Browser snapshot is already filtered to experiments
    that pass these criteria, so all 314 pass by construction. This function is
    therefore a guard for re-derived or newly added data, not a filter that will
    reject anything in the shipped snapshot.
    """
    md = _metadata_public()
    row = md.loc[md["expName"] == expName]
    if row.empty:
        if not has_signed_t(expName):
            raise KeyError(
                f"{expName} is a paper-only experiment (set100/set101); the paper "
                "xlsx does not ship QC fields, so it cannot be QC'd."
            )
        raise KeyError(f"expName not found in QC metadata: {expName}")
    r = row.iloc[0]
    results: dict = {}
    overall = True
    for col, (thresh, op, _desc) in QC_METRICS.items():
        if col not in r.index:
            continue
        val = float(r[col]) if pd.notna(r[col]) else float("nan")
        if op == "ge":
            ok = val >= thresh
        elif op == "le":
            ok = val <= thresh
        else:  # abs_le
            ok = abs(val) <= thresh
        results[col] = (val, thresh, ok)
        overall &= ok
    return overall, results


# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------

def plot_gene_profile(locus: str, out: str | Path, use_paper: bool = True) -> Path:
    """Bar chart of one gene's fitness across all experiments, colored by group.

    Saves a PNG to `out` and returns the Path.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fit = fitness(locus, use_paper=use_paper)
    md = _metadata(use_paper)[["expName", "expGroup", "expDesc"]]
    df = fit.reset_index().merge(md, on="expName", how="left").sort_values(
        ["expGroup", "fit"]
    )

    fig, ax = plt.subplots(figsize=(14, 4))
    groups = df["expGroup"].fillna("(none)").unique()
    palette = plt.cm.tab20.colors  # type: ignore[attr-defined]
    color_map = {g: palette[i % len(palette)] for i, g in enumerate(groups)}
    colors = df["expGroup"].fillna("(none)").map(color_map).tolist()
    ax.bar(range(len(df)), df["fit"], color=colors, width=1.0)
    ax.axhline(0, color="k", linewidth=0.5)
    ax.axhline(-1, color="r", linewidth=0.5, linestyle="--")
    ax.axhline(1, color="r", linewidth=0.5, linestyle="--")
    ax.set_xticks([])
    ax.set_xlabel(f"{len(df)} experiments, grouped by expGroup")
    ax.set_ylabel("fitness (log2 end/start)")
    ax.set_title(f"{locus}  —  fitness profile")
    handles = [plt.Rectangle((0, 0), 1, 1, color=color_map[g]) for g in groups]
    ax.legend(handles, groups, fontsize=7, ncol=4, loc="upper right")
    fig.tight_layout()
    out = Path(out)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def plot_substrate_heatmap(
    grid: pd.DataFrame,
    out: str | Path,
    annot_desc: bool = True,
    vmin: float = -3,
    vmax: float = 3,
) -> Path:
    """Heatmap of a ``substrate_grid()`` result.

    `grid` should be the output of substrate_grid(): index = locusId, experiment
    columns named "<expName> <expDesc>", plus META_COLS carried through.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    cols = [c for c in grid.columns if c not in META_COLS]
    mat = grid[cols].astype(float).to_numpy()

    fig_h = max(4, 0.25 * grid.shape[0])
    fig_w = max(6, 0.6 * len(cols))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    im = ax.imshow(mat, aspect="auto", cmap="RdBu_r", vmin=vmin, vmax=vmax)

    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=45, ha="right", fontsize=7)
    if annot_desc and "sysName" in grid.columns:
        labels = grid["sysName"].astype(str) + "  " + grid["desc"].astype(str).str[:40]
    else:
        labels = grid.index.astype(str)
    ax.set_yticks(range(grid.shape[0]))
    ax.set_yticklabels(labels, fontsize=7)

    cbar = fig.colorbar(im, ax=ax, shrink=0.6)
    cbar.set_label("fitness (log2)")
    fig.tight_layout()
    out = Path(out)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


# ---------------------------------------------------------------------------
# smoke test
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    print(f"data_dir: {data_dir()}")
    print("\nSmoke test on PP_0154 (propionyl-CoA:succinate CoA transferase)")
    g = gene("PP_0154")
    print(f"  desc: {g.get('desc')}")
    print(f"  mean_fit: {g['mean_fit']:.3f}, min: {g['min_fit']:.3f}, max: {g['max_fit']:.3f}")
    print(f"  n_phenotypes: {g['n_phenotypes']}")

    print("\nTop 5 phenotypes:")
    print(phenotypes("PP_0154").head().to_string(index=False))

    print("\nTop 5 cofit neighbors:")
    print(cofitness("PP_0154", n=5).to_string(index=False))

    print("\nAcetate experiments (paper):")
    ace = experiments(contains="acetate")
    print(ace[["expName", "expDesc", "expGroup"]].head().to_string(index=False))
