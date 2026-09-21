#!/usr/bin/env python
"""
substrate_map.py — heatmap of genes × substrates.

Usage:
    python substrate_map.py PP_0154,PP_4487,PP_1743 acetate,glucose,citrate
    python substrate_map.py PP_0154,PP_4487 --group "carbon source" --out /tmp/cmap.png
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

SCRIPTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SCRIPTS))

from putida_rbtnseq import (  # noqa: E402
    substrate_grid,
    plot_substrate_heatmap,
    experiments,
)


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("genes", help="comma-separated locusIds, e.g. PP_0154,PP_4487")
    p.add_argument(
        "substrates",
        nargs="?",
        default=None,
        help="comma-separated substrate substrings (or omit and use --group)",
    )
    p.add_argument("--group", help="filter experiments by expGroup (e.g. 'carbon source')")
    p.add_argument("--out", default="/tmp/substrate_map.png", help="output PNG path")
    args = p.parse_args()

    genes = [g.strip() for g in args.genes.split(",") if g.strip()]
    if args.substrates:
        subs = [s.strip() for s in args.substrates.split(",") if s.strip()]
    elif args.group:
        md = experiments(group=args.group)
        # use unique expDesc values as filter substrings
        subs = md["expDesc"].dropna().unique().tolist()
        print(f"Using {len(subs)} experiments from group '{args.group}'")
    else:
        sys.exit("Provide substrates (positional) or --group")

    grid = substrate_grid(genes=genes, substrates=subs)
    print(f"Grid shape: {grid.shape[0]} genes × {len([c for c in grid.columns if c not in ['orgId','sysName','geneName','desc']])} experiments")

    # quick text summary before plotting
    fitcols = [c for c in grid.columns if c not in ("orgId", "sysName", "geneName", "desc")]
    print("\n--- fitness values ---")
    print(grid[["sysName", "desc"] + fitcols[: min(8, len(fitcols))]].to_string())

    out = plot_substrate_heatmap(grid, out=args.out)
    print(f"\nHeatmap → {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
