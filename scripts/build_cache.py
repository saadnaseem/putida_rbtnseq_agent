"""
Build the parquet cache for the P. putida KT2440 RB-TnSeq skill.

Reads the raw TSV/XLSX/CSV source files and writes parquet to ``<data>/cache/``.
The helper module (``putida_rbtnseq.py``) reads only the cache, never the originals.

Data location is resolved as:
    1. $PUTIDA_RBTNSEQ_DATA
    2. <repo root>/data

See docs/DATA.md for how to obtain the source files.

Usage:
    python scripts/build_cache.py
    python scripts/build_cache.py --rebuild
    python scripts/build_cache.py --manifest      # write cache/MANIFEST.json only
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from putida_rbtnseq import cache_dir, data_dir  # noqa: E402

META_COLS = ["orgId", "locusId", "sysName", "geneName", "desc"]

# Source files. Each entry maps a logical name to the accepted filename patterns —
# browser downloads often append " (1)", so more than one spelling is allowed.
SOURCES: dict[str, list[str]] = {
    "fit_tsv":      ["fit_organism_Putida.tsv", "fit_organism_Putida (*).tsv"],
    "t_tsv":        ["t_organism_Putida.tsv", "t_organism_Putida (*).tsv"],
    "exp_txt":      ["exp_organism_Putida.txt"],
    "genes_tab":    ["organism_Putida_genes.tab"],
    "reanno_tsv":   ["reanno_Putida.tsv"],
    "cofit_txt":    ["cofit_organism_Putida.txt"],
    "specific_txt": ["specific_phenotypes_Putida.txt"],
    "uniprot_csv":  ["Pseudomonas_putida_protein_list_*.csv", "uniprot_putida.csv"],
    "fmodule_xlsx": ["fModule_Metadata.xlsx"],
}


def resolve(key: str) -> Path:
    """Find a source file by any of its accepted names."""
    d = data_dir()
    for pattern in SOURCES[key]:
        hits = sorted(d.glob(pattern))
        if hits:
            return hits[0]
    raise SystemExit(
        f"Missing source file for '{key}' in {d}\n"
        f"  expected one of: {SOURCES[key]}\n"
        f"  see docs/DATA.md for download instructions"
    )


def normalize_object_dtypes(df: pd.DataFrame) -> pd.DataFrame:
    """Cast object columns to string so pyarrow serializes them cleanly.
    The xlsx export has columns like `dateStarted` with mixed datetime/str."""
    for col in df.select_dtypes(include=["object", "string"]).columns:
        df[col] = df[col].astype("string")
    return df


def melt_matrix(df: pd.DataFrame, value_name: str) -> pd.DataFrame:
    """Wide gene x experiment matrix -> long (locusId, expName, value)."""
    long = df.melt(id_vars=META_COLS, var_name="expCol", value_name=value_name)
    # column headers look like "set10IT006 L-Arginine (C)" — split on the first space
    long["expName"] = long["expCol"].str.split(" ", n=1).str[0]
    long["expDesc"] = long["expCol"].str.split(" ", n=1).str[1]
    return long[["locusId", "sysName", "expName", "expDesc", value_name]]


def _melt_and_join(fit: pd.DataFrame, t: pd.DataFrame) -> pd.DataFrame:
    fit_long = melt_matrix(fit, "fit")
    t_long = melt_matrix(t, "t")[["locusId", "expName", "t"]]
    long = fit_long.merge(t_long, on=["locusId", "expName"], how="left")
    expected = fit.shape[0] * (fit.shape[1] - len(META_COLS))
    if len(long) != expected:
        raise SystemExit(
            f"melt/join produced {len(long):,} rows, expected {expected:,} — "
            "duplicate expName columns in the source matrices?"
        )
    return normalize_object_dtypes(long)


def build_fitness(rebuild: bool) -> None:
    C = cache_dir()
    outs = [C / "fit.parquet", C / "t.parquet", C / "fit_long.parquet"]
    if not rebuild and all(p.exists() for p in outs):
        print("[skip] public fit/t parquet already present")
        return

    fit_tsv, t_tsv = resolve("fit_tsv"), resolve("t_tsv")
    print(f"[load] {fit_tsv.name}")
    fit = normalize_object_dtypes(pd.read_csv(fit_tsv, sep="\t"))
    print(f"       {fit.shape[0]} genes x {fit.shape[1] - len(META_COLS)} experiments")
    fit.to_parquet(outs[0], index=False)

    print(f"[load] {t_tsv.name}")
    t = normalize_object_dtypes(pd.read_csv(t_tsv, sep="\t"))
    t.to_parquet(outs[1], index=False)

    print("[melt] public fit + t -> long")
    long = _melt_and_join(fit, t)
    long.to_parquet(outs[2], index=False)
    print(f"       {len(long):,} rows")


def build_xlsx_matrices(rebuild: bool) -> None:
    """Same matrices, but the paper version (332 experiments) — canonical."""
    C = cache_dir()
    outs = [
        C / "fit_paper.parquet", C / "t_paper.parquet",
        C / "fit_paper_long.parquet", C / "metadata_paper.parquet",
    ]
    if not rebuild and all(p.exists() for p in outs):
        print("[skip] paper xlsx parquet already present")
        return

    xlsx = resolve("fmodule_xlsx")
    print(f"[load] {xlsx.name} — three sheets, this is slow (~30 s)")
    fit = normalize_object_dtypes(pd.read_excel(xlsx, sheet_name="fitness_measurements"))
    t = normalize_object_dtypes(pd.read_excel(xlsx, sheet_name="T-like_statistics"))
    meta = normalize_object_dtypes(pd.read_excel(xlsx, sheet_name="metadata"))
    print(f"       fit  {fit.shape}\n       t    {t.shape}\n       meta {meta.shape}")

    fit.to_parquet(outs[0], index=False)
    t.to_parquet(outs[1], index=False)
    meta.to_parquet(outs[3], index=False)

    print("[melt] paper fit + t -> long")
    long = _melt_and_join(fit, t)
    long.to_parquet(outs[2], index=False)
    print(f"       {len(long):,} rows")


def _simple(key: str, out_name: str, rebuild: bool, sep: str | None, label: str,
            rename: dict[str, str] | None = None) -> None:
    pq = cache_dir() / out_name
    if not rebuild and pq.exists():
        print(f"[skip] {out_name} already present")
        return
    src = resolve(key)
    print(f"[load] {src.name}")
    df = pd.read_csv(src, sep=sep) if sep else pd.read_csv(src)
    if rename:
        df = df.rename(columns=rename)
    df = normalize_object_dtypes(df)
    df.to_parquet(pq, index=False)
    print(f"       {len(df):,} {label}")


def write_manifest() -> None:
    """Record row/column counts and checksums so a rebuild can be verified."""
    C = cache_dir()
    entries = {}
    for f in sorted(C.glob("*.parquet")):
        df = pd.read_parquet(f)
        h = hashlib.sha256(f.read_bytes()).hexdigest()[:16]
        entries[f.name] = {
            "rows": int(df.shape[0]),
            "cols": int(df.shape[1]),
            "bytes": f.stat().st_size,
            "sha256_16": h,
        }
    manifest = {
        "built_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "pandas": pd.__version__,
        "files": entries,
    }
    (C / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\n[manifest] {C / 'MANIFEST.json'}")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--rebuild", action="store_true", help="force rebuild of all caches")
    p.add_argument("--manifest", action="store_true", help="only rewrite MANIFEST.json")
    args = p.parse_args()

    D = data_dir()
    if not D.is_dir():
        sys.exit(
            f"Data directory not found: {D}\n"
            "Set $PUTIDA_RBTNSEQ_DATA or create <repo>/data — see docs/DATA.md"
        )
    cache_dir().mkdir(exist_ok=True, parents=True)
    print(f"data_dir : {D}\ncache_dir: {cache_dir()}\n")

    if args.manifest:
        write_manifest()
        return 0

    _simple("genes_tab", "genes.parquet", args.rebuild, "\t", "genes")
    _simple("reanno_tsv", "reanno.parquet", args.rebuild, "\t", "curated reannotations")
    _simple("uniprot_csv", "uniprot.parquet", args.rebuild, None, "UniProt rows",
            rename={"Gene Names (ordered locus)": "locusId"})
    _simple("exp_txt", "metadata.parquet", args.rebuild, "\t", "experiments")
    build_fitness(args.rebuild)
    build_xlsx_matrices(args.rebuild)
    _simple("cofit_txt", "cofit.parquet", args.rebuild, "\t", "cofitness pairs")
    _simple("specific_txt", "specific_phenotypes.parquet", args.rebuild, "\t",
            "significant gene-condition calls")

    print(f"\nCache built at {cache_dir()}")
    for f in sorted(cache_dir().glob("*.parquet")):
        print(f"  {f.name:35s}  {f.stat().st_size / 1e6:6.2f} MB")
    write_manifest()
    return 0


if __name__ == "__main__":
    sys.exit(main())
