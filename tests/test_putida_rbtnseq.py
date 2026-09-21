"""Regression tests for the P. putida RB-TnSeq helper API.

Run:  pytest -q

Tests requiring the dataset skip automatically when the cache is absent, so the
suite is safe to run in CI without redistributing the data.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import putida_rbtnseq as P  # noqa: E402

# --------------------------------------------------------------------------

def _cache_present() -> bool:
    return (P.cache_dir() / "fit_paper_long.parquet").exists()


needs_data = pytest.mark.skipif(
    not _cache_present(),
    reason=f"parquet cache not found at {P.cache_dir()} — see docs/DATA.md",
)

# Reference gene: succinyl-CoA:acetate CoA-transferase, acetate catabolism
GENE = "PP_0154"


# ---------------------------------------------------------------- data shape

@needs_data
def test_matrix_dimensions():
    """Paper matrix = 332 experiments, public snapshot = 314."""
    assert P._long_paper()["expName"].nunique() == 332
    assert P._long_public()["expName"].nunique() == 314
    assert P._fit_paper().shape == (4732, 337)   # 5 meta + 332 experiments
    assert P._fit_public().shape == (4778, 319)  # 5 meta + 314 experiments


@needs_data
def test_no_duplicate_gene_experiment_pairs():
    for long in (P._long_paper(), P._long_public()):
        assert long.duplicated(["locusId", "expName"]).sum() == 0


@needs_data
def test_metadata_covers_every_experiment():
    assert set(P._long_paper()["expName"]) <= set(P._metadata_paper()["expName"])
    assert set(P._long_public()["expName"]) <= set(P._metadata_public()["expName"])


@needs_data
def test_no_missing_fitness_values():
    long = P._long_paper()
    assert long["fit"].isna().sum() == 0
    assert long["t"].isna().sum() == 0


# ------------------------------------------------- the unsigned-t data caveat

@needs_data
def test_set100_and_set101_store_unsigned_t():
    """Regression guard for the documented unsigned-t caveat (docs/DATA_CAVEATS.md).

    In the paper xlsx, set100 and set101 report `t` as a magnitude: zero negative
    values. Every other set is roughly half negative, as a signed statistic must be.
    """
    long = P._long_paper().assign(
        setid=lambda d: d["expName"].str.extract(r"^(set\d+)")[0]
    )
    for setid, grp in long.groupby("setid"):
        frac_neg = (grp["t"] < 0).mean()
        if setid in P.SETS_WITHOUT_SIGNED_T:
            assert (grp["t"] < 0).sum() == 0, f"{setid} unexpectedly has signed t"
        else:
            assert 0.4 < frac_neg < 0.6, f"{setid} has {frac_neg:.2f} negative t"


@needs_data
def test_public_snapshot_is_sign_consistent():
    """The public TSV has no sign anomaly: fit and t always agree in sign."""
    long = P._long_public().dropna(subset=["fit", "t"])
    sig = long[(long["fit"].abs() > 1) & (long["t"].abs() > 4)]
    assert (np.sign(sig["fit"]) != np.sign(sig["t"])).sum() == 0


def test_has_signed_t_flag():
    assert P.has_signed_t("set1IT084")
    assert not P.has_signed_t("set100IT039")
    assert not P.has_signed_t("set101IT012")


# ------------------------------------------------------- end-to-end validation

@needs_data
def test_significance_rule_recovers_fitness_browser_calls():
    """|fit|>1 & |t|>4 on the public matrix recovers the Fitness Browser's own
    specific-phenotype calls. This validates the whole load -> melt -> join path."""
    spec = P._specific()
    fb = set(zip(spec["locusId"].astype(str), spec["expName"].astype(str)))
    pub = P._long_public().dropna(subset=["fit", "t"])
    mine = pub[(pub["fit"].abs() > 1) & (pub["t"].abs() > 4)]
    ours = set(zip(mine["locusId"].astype(str), mine["expName"].astype(str)))
    recall = len(fb & ours) / len(fb)
    assert recall > 0.999, f"recall of Fitness Browser calls dropped to {recall:.3f}"


# -------------------------------------------------------------- gene-level API

@needs_data
def test_gene_happy_path():
    g = P.gene(GENE)
    assert g["has_fitness_data"] is True
    assert g["n_experiments"] == 332
    assert "CoA transferase" in str(g["desc"])
    assert g["n_phenotypes"] > 0


@needs_data
def test_gene_without_fitness_data_does_not_raise():
    """929 of 5,661 annotated genes have no fitness data (essential / no insertions).
    gene() must return annotation for them rather than raising."""
    scored = set(P._long_paper()["locusId"])
    unscored = sorted(set(P._genes()["locusId"]) - scored)
    assert len(unscored) > 0
    g = P.gene(unscored[0])
    assert g["has_fitness_data"] is False
    assert g["n_experiments"] == 0
    assert g["n_phenotypes"] == 0
    assert pd.isna(g["mean_fit"])
    assert str(g["desc"])  # annotation still present


@needs_data
def test_unknown_gene_raises():
    with pytest.raises(KeyError):
        P.gene("PP_99999")


@needs_data
def test_phenotype_columns_are_stable_when_empty():
    """Empty and non-empty results must share column order so they concatenate."""
    non_empty = P.phenotypes(GENE)
    assert list(non_empty.columns) == P.PHENOTYPE_COLS
    assert len(non_empty) > 0

    scored = P._long_paper()["locusId"].unique()
    empty = next(
        (P.phenotypes(l) for l in scored[:500] if len(P.phenotypes(l)) == 0), None
    )
    assert empty is not None, "expected at least one gene with no phenotypes"
    assert list(empty.columns) == P.PHENOTYPE_COLS
    pd.concat([non_empty, empty])  # must not raise or misalign


@needs_data
def test_phenotypes_respect_thresholds():
    ph = P.phenotypes(GENE)
    assert (ph["fit"].abs() > P.FIT_THRESH).all()
    assert (ph["t"].abs() > P.T_THRESH).all()
    assert ph["fit"].is_monotonic_increasing  # sorted most-negative first


@needs_data
def test_known_acetate_phenotype():
    """Biological anchor: PP_0154 is specifically required on acetate."""
    ph = P.phenotypes(GENE)
    row = ph[ph["expName"] == "set1IT084"]
    assert len(row) == 1
    assert row["fit"].iloc[0] == pytest.approx(-1.521, abs=0.01)
    assert abs(row["t"].iloc[0]) > 9


# ------------------------------------------------------------------ cofitness

@needs_data
def test_cofitness_sorted_by_rank():
    co = P.cofitness(GENE, n=10)
    assert co["rank"].is_monotonic_increasing
    assert co["cofit"].is_monotonic_decreasing
    assert co["hitId"].iloc[0] == "PP_4487"  # acetyl-CoA synthetase


# ------------------------------------------------------------------------- QC

@needs_data
def test_every_public_experiment_passes_qc():
    """The public snapshot is pre-filtered by the Fitness Browser to experiments
    meeting the Wetmore 2015 criteria, so all 314 must pass."""
    md = P._metadata_public()
    failures = [e for e in md["expName"] if not P.qc_pass(e)[0]]
    assert failures == [], f"{len(failures)} experiments failed QC: {failures[:5]}"


@needs_data
def test_qc_metrics_use_documented_columns():
    """Guard against the mad12/mad12c mix-up: mad12c is a correlation-like value
    centred near 0.96 and would reject every experiment under a <=0.5 threshold."""
    assert "mad12" in P.QC_METRICS
    assert "mad12c" not in P.QC_METRICS
    md = P._metadata_public()
    assert md["mad12"].max() <= 0.5
    assert md["mad12c"].min() > 0.5  # the column that must NOT be used


@needs_data
def test_qc_on_paper_only_experiment_explains_itself():
    with pytest.raises(KeyError, match="paper-only"):
        P.qc_pass("set100IT039")


# ----------------------------------------------------------------- conditions

@needs_data
def test_experiments_filters():
    assert len(P.experiments(contains="acetate")) == 4
    assert len(P.experiments(group="carbon source")) > 100
    assert len(P.experiments(contains="zzzznotathing")) == 0


@needs_data
def test_condition_top_genes_signs():
    neg = P.condition_top_genes("set1IT084", n=10, sign="neg")
    assert (neg["fit"] < -P.FIT_THRESH).all()
    pos = P.condition_top_genes("set1IT084", n=10, sign="pos")
    assert (pos["fit"] > P.FIT_THRESH).all()
    with pytest.raises(ValueError):
        P.condition_top_genes("set1IT084", sign="sideways")
    with pytest.raises(KeyError):
        P.condition_top_genes("not_an_experiment")


@needs_data
def test_condition_compare():
    cmp = P.condition_compare("set1IT084", "set1IT078")
    assert list(cmp.columns) == [
        "locusId", "sysName", "desc", "fit_a", "t_a", "fit_b", "t_b", "delta"
    ]
    assert len(cmp) > 0
    assert cmp["delta"].abs().is_monotonic_decreasing
    np.testing.assert_allclose(cmp["delta"], cmp["fit_a"] - cmp["fit_b"], atol=1e-9)


# -------------------------------------------------------------------- grid

@needs_data
def test_substrate_grid():
    grid = P.substrate_grid([GENE, "PP_4487", "PP_1743"], ["acetate"])
    assert grid.index.name == "locusId"
    assert len(grid) == 3
    exp_cols = [c for c in grid.columns if c not in P.META_COLS]
    assert all("acetate" in c.lower() for c in exp_cols)
    with pytest.raises(ValueError):
        P.substrate_grid([GENE], ["zzzznotathing"])


# ------------------------------------------------------------------- plotting

@needs_data
def test_plots_write_files(tmp_path):
    p1 = P.plot_gene_profile(GENE, tmp_path / "profile.png")
    assert p1.exists() and p1.stat().st_size > 1000
    grid = P.substrate_grid([GENE, "PP_4487"], ["acetate"])
    p2 = P.plot_substrate_heatmap(grid, tmp_path / "heat.png")
    assert p2.exists() and p2.stat().st_size > 1000


# --------------------------------------------------------- data-absent behaviour

def test_missing_cache_gives_actionable_error(monkeypatch, tmp_path):
    monkeypatch.setenv("PUTIDA_RBTNSEQ_DATA", str(tmp_path))
    P._genes.cache_clear()
    try:
        with pytest.raises(P.DataNotFound, match="docs/DATA.md"):
            P._genes()
    finally:
        P._genes.cache_clear()
