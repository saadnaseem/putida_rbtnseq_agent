# fModules — context for the Borchert 2024 paper

## What an fModule is

Borchert, Bleem, Lim et al. 2024 (mSystems) ran **Independent Component Analysis (ICA)** on a 4,732 × 332 fitness matrix for *P. putida* KT2440 (the same `fitness_measurements` sheet in `fModule_Metadata.xlsx`). ICA decomposes the matrix into:

```
fitness  ≈  M  ×  A
(genes×expts)  (genes×k)  (k×expts)
```

where `k=84` is the chosen number of independent components. Each component is interpreted as an **fModule** = a coherent group of genes that vary together in fitness across experiments. Genes with strong loading (|weight| > some threshold) on a component are "members" of that fModule.

The paper validated 4 modules with engineered mutants — hydroxycinnamate metabolism, acetyl-CoA assimilation, nitrogen metabolism, and a stress-resistance module.

## What we have vs what we don't

We **have**:
- The input fitness matrix (`fit_paper.parquet`)
- The input t-stat matrix (`t_paper.parquet`)
- The experiment metadata (`metadata_paper.parquet`)
- The paper PDF

We **don't have**:
- The 84 ICA mixing matrix (`M`) — i.e. the gene-by-module loadings
- The fModule labels / interpretations

## How to use the concept anyway

### Option 1 — `cofitness()` as a pairwise proxy

For most "what's in the same module as gene X" questions, the precomputed cofitness table gives ~90% of the answer:

```python
from putida_rbtnseq import cofitness
print(cofitness("PP_0154", n=20))
```

A cohesive set of genes with high mutual cofitness ≈ an fModule.

### Option 2 — Re-derive ICA modules locally (~30 s)

```python
import sys
sys.path.insert(0, "$SKILL/scripts")
from putida_rbtnseq import _fit_paper, META_COLS
import pandas as pd
from sklearn.decomposition import FastICA

fit = _fit_paper().set_index("locusId")
mat = fit.drop(columns=[c for c in META_COLS if c in fit.columns]).fillna(0).to_numpy()

ica = FastICA(n_components=84, random_state=0, max_iter=500, tol=1e-4)
ica.fit(mat)
M = pd.DataFrame(ica.components_.T, index=fit.index)   # genes × 84 modules
print(M.shape)

# Find modules where PP_0154 has highest absolute loading
import numpy as np
loadings = M.loc["PP_0154"]
top_modules = loadings.abs().sort_values(ascending=False).head(5).index
for m in top_modules:
    members = M[m].abs().sort_values(ascending=False).head(15).index.tolist()
    print(f"\nModule {m}, PP_0154 loading {loadings[m]:.2f}")
    print("  members:", ", ".join(members))
```

Note: ICA is sign-arbitrary and seed-sensitive — your modules won't have the same numeric IDs as the paper's, but the gene groupings should largely match. For a stable run, set `random_state` and use the paper's `n_components=84`.

### Option 3 — Live correlation against the full matrix

For a single query gene's "module-mates":

```python
import sys
sys.path.insert(0, "$SKILL/scripts")
from putida_rbtnseq import _fit_paper, META_COLS

fit = _fit_paper().set_index("locusId")
mat = fit.drop(columns=[c for c in META_COLS if c in fit.columns]).fillna(0)

target = mat.loc["PP_0154"]
corr = mat.apply(lambda row: row.corr(target), axis=1)
print(corr.sort_values(ascending=False).head(20))
```

This is essentially what `cofit_organism_Putida.txt` was generated from.

## When fModules matter

- **Hypothesis generation**: "PP_X is in a module with established acetate-uptake genes → likely involved in acetate metabolism"
- **Functional annotation**: an unannotated gene grouped with biosynthesis genes for amino acid Y is probably also involved in Y biosynthesis
- **Pathway validation**: do all genes you'd expect in a pathway cluster together? If not, why not? Polar effects? Redundancy? Multiple isoforms?

If the user asks specifically about *the paper's* modules by number (e.g. "module 14"), tell them we don't have the loadings table — point them to the paper's supplementary materials.
