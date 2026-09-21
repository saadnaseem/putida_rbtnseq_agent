#!/usr/bin/env python
"""
condition_phenotypes.py — top phenotypes for a single condition.

Usage:
    python condition_phenotypes.py set1IT084
    python condition_phenotypes.py set1IT084 --sign pos
    python condition_phenotypes.py set1IT084 --top 30 --sign abs

If you don't know the expName, run:
    python -c "import sys; sys.path.insert(0,'..'); from putida_rbtnseq import experiments; \
        print(experiments(contains='acetate')[['expName','expDesc']].to_string(index=False))"
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))

from putida_rbtnseq import (  # noqa: E402
    condition_top_genes,
    experiments,
    qc_pass,
    _metadata_paper,
)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("expName", help="experiment ID, e.g. set1IT084")
    p.add_argument("--top", type=int, default=20)
    p.add_argument("--sign", choices=["neg", "pos", "abs"], default="neg",
                   help="neg=depleted (required), pos=enriched (inhibitory), abs=either")
    args = p.parse_args()

    md = _metadata_paper()
    row = md.loc[md["expName"] == args.expName]
    if row.empty:
        print(f"expName not found: {args.expName}")
        sims = experiments(contains=args.expName)
        if not sims.empty:
            print("Did you mean one of:")
            print(sims[["expName", "expDesc"]].head(5).to_string(index=False))
        return 1

    r = row.iloc[0]
    print(f"=== {args.expName} ===")
    print(f"  desc:  {r.get('expDesc')}")
    print(f"  group: {r.get('expGroup')}")
    print(f"  media: {r.get('media')}   condition: {r.get('condition_1')} {r.get('units_1')} {r.get('concentration_1')}")

    try:
        ok, _ = qc_pass(args.expName)
        print(f"  QC:    {'PASS' if ok else 'FAIL'}  (call qc_pass() for details)")
    except KeyError:
        print(f"  QC:    (not in public-snapshot metadata — can't assess)")

    print(f"\n--- Top {args.top} genes (sign={args.sign}) ---")
    df = condition_top_genes(args.expName, n=args.top, sign=args.sign)
    if df.empty:
        print("  (no significant genes at default thresholds)")
    else:
        print(df.to_string(index=False))

    return 0


if __name__ == "__main__":
    sys.exit(main())
