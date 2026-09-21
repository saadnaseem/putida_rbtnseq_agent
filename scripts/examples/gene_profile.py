#!/usr/bin/env python
"""
gene_profile.py — one-page fitness profile for a P. putida gene.

Usage:
    python gene_profile.py PP_0154
    python gene_profile.py PP_0154 --plot /tmp/PP_0154.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))

import pandas as pd  # noqa: E402

from putida_rbtnseq import (  # noqa: E402
    gene,
    phenotypes,
    cofitness,
    plot_gene_profile,
)


def _val(d, key):
    """Return d[key] if it's a non-null, non-empty value, else None."""
    v = d.get(key)
    if v is None or (isinstance(v, float) and v != v):  # NaN
        return None
    try:
        if pd.isna(v):
            return None
    except (TypeError, ValueError):
        pass
    s = str(v).strip()
    return s if s and s.lower() not in ("nan", "<na>") else None


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("locus", help="locusId, e.g. PP_0154")
    p.add_argument("--top", type=int, default=10, help="top-N phenotypes / cofit hits")
    p.add_argument("--plot", help="save profile PNG to this path")
    args = p.parse_args()

    g = gene(args.locus)
    title = _val(g, "reannotation") or _val(g, "Protein names") or _val(g, "desc") or "?"

    print(f"=== {args.locus}  —  {title} ===")
    ec = _val(g, "EC number")
    if ec:
        print(f"  EC: {ec}")
    pw = _val(g, "Pathway")
    if pw:
        print(f"  Pathway: {pw[:120]}")
    print(f"  scaffold: {_val(g, 'scaffoldId')}  pos: {g.get('begin')}–{g.get('end')} ({_val(g, 'strand')})")
    print(f"  fitness:  mean={g['mean_fit']:.2f}  min={g['min_fit']:.2f}  max={g['max_fit']:.2f}")
    print(f"  n_phenotypes (|fit|>1, |t|>4): {g['n_phenotypes']}")

    print("\n--- Top phenotypes (signed by fit) ---")
    ph = phenotypes(args.locus)
    if ph.empty:
        print("  (no significant phenotypes)")
    else:
        # collapse replicates by expDesc
        collapsed = (
            ph.groupby(["expGroup", "expDesc"], as_index=False)
              .agg(fit=("fit", "mean"),
                   t_abs_max=("t", lambda s: float(s.abs().max())),
                   n_replicates=("expName", "size"))
              .sort_values("fit")
        )
        print(collapsed.head(args.top).to_string(index=False))

    print("\n--- Top cofit neighbors ---")
    co = cofitness(args.locus, n=args.top)
    print(co.to_string(index=False))

    if args.plot:
        out = plot_gene_profile(args.locus, args.plot)
        print(f"\nProfile chart → {out}")

    # 3-line summary
    print("\n--- Summary ---")
    print(f"{args.locus} — {title}")
    if not ph.empty:
        top = ph.iloc[0]
        print(f"  · Strongest phenotype: {top['expDesc']} (fit={top['fit']:.2f}, |t|={abs(top['t']):.1f})")
    if not co.empty:
        nb = co.iloc[0]
        hd = str(nb['hitDesc'])[:60] if pd.notna(nb['hitDesc']) else nb['hitSysName']
        print(f"  · Top cofit: {nb['hitId']} {hd} (cofit={nb['cofit']:.2f})")
    print(f"  · Browser: https://fit.genomics.lbl.gov/cgi-bin/myFitShow.cgi?orgId=Putida&gene={args.locus}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
