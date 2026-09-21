---
name: putida-rbtnseq
description: Query and analyze RB-TnSeq gene-fitness data for *Pseudomonas putida* KT2440 from the LBNL Fitness Browser + Borchert 2024 fModule paper. Use when the user asks about a P. putida gene (locusId starting with `PP_`), a substrate-utilization phenotype, a fitness profile across conditions, cofitness/co-essentiality, fModules, or compares conditions (e.g. acetate vs glucose). Triggers on: "P. putida", "Pseudomonas putida", "KT2440", "RB-TnSeq", "fitness browser", "fitness data", "PP_NNNN", "fModule", "cofit", "specific phenotype", "substrate utilization", "carbon source phenotype". Do NOT trigger for other organisms, RNA-seq/proteomics analysis, or DBTL plate-reader cleanup.
---

# P. putida RB-TnSeq skill

Random-Barcoded Transposon-Sequencing fitness data for *P. putida* KT2440. Canonical
source: paper version (332 experiments, Borchert 2024). Public Fitness Browser snapshot
(314 experiments) also available via `use_paper=False`.

## Setup

The dataset is **not** bundled. If the cache is missing, every call raises `DataNotFound`
with instructions. To set it up, follow `docs/DATA.md` (download ~120 MB, then
`python scripts/build_cache.py`).

```bash
SKILL=$(dirname "$0")                      # this skill's directory
export PUTIDA_RBTNSEQ_DATA=~/data/rbtnseq  # where the dataset lives
```

## How to invoke

All analysis goes through one helper module — **never read the raw TSV/XLSX**, the
parquet cache loads in <50 ms.

Either run a one-liner:

```bash
python -c "import sys; sys.path.insert(0, '$SKILL/scripts'); from putida_rbtnseq import phenotypes; print(phenotypes('PP_0154').to_string(index=False))"
```

or run an example script:

```bash
python $SKILL/scripts/examples/gene_profile.py PP_0154
python $SKILL/scripts/examples/condition_phenotypes.py set1IT084
python $SKILL/scripts/examples/substrate_map.py PP_0154,PP_4487,PP_1743 acetate,glucose,citrate
```

## Decision tree

| User asks about… | Call this | Plus |
|---|---|---|
| One specific gene (`PP_NNNN`) | `gene(locus)` → `phenotypes(locus)` → `cofitness(locus)` | `plot_gene_profile(locus, out=...)` for a chart |
| What gene does in one condition | `condition_top_genes(expName, sign='neg')` | `experiments(contains='glucose')` to find expNames first |
| Which conditions test substrate X | `experiments(contains='X')` | then loop `condition_top_genes` over them |
| Substrate-utilization map (genes × substrates) | `substrate_grid(genes, substrates)` | `plot_substrate_heatmap(grid, out=...)` |
| Compare two conditions | `condition_compare(exp_a, exp_b)` | sort by `delta` for biggest swings |
| Pathway / module of related genes | `cofitness(locus, n=50)` | see `references/fmodules.md` for the ICA paper context |
| Is this experiment trustworthy? | `qc_pass(expName)` | public-snapshot experiments only |

## Standard significance call

A gene–condition phenotype is "specific" when **|fit| > 1 AND |t| > 4**. Negative `fit`
= mutant depleted = gene important for growth. Helper functions use these defaults;
override with `fit_thresh=` / `t_thresh=`.

## Two things that will produce wrong answers

1. **Never interpret `sign(t)`.** `set100` and `set101` (78 of 332 experiments) store `t`
   as an unsigned magnitude — a negative `fit` with a positive `t` is expected there, not
   a contradiction. Use signed `fit` for direction, `|t|` for confidence. Check with
   `has_signed_t(expName)`.
2. **929 of 5,661 genes have no fitness data** (essential, or too few insertions).
   `gene()` returns `has_fitness_data=False` for them. Absence of a phenotype is a
   measurement gap, not evidence the gene is dispensable — say so rather than reporting
   "no phenotype".

Full detail in `docs/DATA_CAVEATS.md`.

## When to read references (lazy load only what's needed)

- `references/files.md` — what each raw data file contains, schema, join keys
- `references/biology.md` — RB-TnSeq primer, fitness/t interpretation, QC thresholds
- `references/workflows.md` — end-to-end recipes (gene profile, substrate panel, condition compare, pathway sanity check)
- `references/fmodules.md` — Borchert 2024 ICA / fModule context
- `docs/DATA_CAVEATS.md` — the unsigned-`t` sets and the QC-column trap

## Always end gene/condition analyses with

A 3-line summary in this format:

```
PP_0154 — succinyl-CoA:acetate CoA-transferase (reanno: same)
  · Specifically required for acetate growth (set1IT084: fit=-1.52, |t|=9.9)
  · Top cofit: PP_4487 acetyl-CoA synthetase (cofit=0.57) → acetate-uptake module
```

This makes it easy to scan when answering iteratively. Report `|t|`, not signed `t`.

## Cross-references

- **Browser**: https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=Putida
  (deep-link a gene with `myFitShow.cgi?orgId=Putida&gene=PP_NNNN`)
- **Paper**: Borchert et al. 2024, *mSystems* — https://doi.org/10.1128/msystems.00934-24
- **Data setup**: `docs/DATA.md`
