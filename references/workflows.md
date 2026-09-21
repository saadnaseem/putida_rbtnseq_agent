# Analysis recipes

All snippets assume:

```python
import sys
sys.path.insert(0, "$SKILL/scripts")
from putida_rbtnseq import (
    gene, fitness, phenotypes, cofitness, experiments,
    condition_top_genes, condition_compare, substrate_grid,
    qc_pass, plot_gene_profile, plot_substrate_heatmap,
)
```

Run with `python`.

---

## 1. Gene profile (the most common request)

> "What does PP_0154 do, and where does it matter?"

```python
g = gene("PP_0154")
print(f"{g['locusId']}  {g.get('reannotation', g['desc'])}")

# top phenotypes — collapsing replicates
ph = phenotypes("PP_0154")
top = (
    ph.groupby("expDesc", as_index=False)
      .agg(fit=("fit", "mean"), t_max=("t", lambda s: s.abs().max()), n=("expName", "size"))
      .sort_values("fit")
      .head(10)
)
print(top.to_string(index=False))

# top cofit neighbors
print(cofitness("PP_0154", n=10).to_string(index=False))

# optional plot
plot_gene_profile("PP_0154", out="/tmp/PP_0154_profile.png")
```

End with the 3-line summary format from SKILL.md.

---

## 2. Substrate-utilization map

> "Build a heatmap of central-carbon genes across all carbon-source experiments."

```python
# pick a gene set (here: TCA cycle in P. putida)
genes = ["PP_0154", "PP_4011", "PP_4012", "PP_4013", "PP_1755", "PP_1756"]

# pick experiments — by group
carbon = experiments(group="carbon source")
exp_descs = carbon["expDesc"].unique().tolist()

grid = substrate_grid(genes=genes, substrates=exp_descs)
plot_substrate_heatmap(grid, out="/tmp/tca_carbon_heatmap.png")
```

If the substrate list is short and named, pass them directly:

```python
grid = substrate_grid(
    genes=["PP_0154", "PP_4487", "PP_1743", "PP_1635"],
    substrates=["acetate", "glucose", "citrate", "succinate"],
)
```

`substrate_grid()` matches substrates as case-insensitive substrings against the experiment column headers (`<expName> <expDesc>`).

---

## 3. What's important for growth in condition X

> "Which genes are required to grow on acetate?"

```python
# 1. find the experiments
ace = experiments(contains="acetate")
print(ace[["expName", "expDesc"]].to_string(index=False))

# 2. top depleted (= required) genes per experiment
for exp in ace["expName"].head(2):
    print(f"\n=== {exp} ===")
    print(condition_top_genes(exp, n=10, sign="neg").to_string(index=False))
```

For a more robust list, intersect across replicates: take genes that score significantly in ≥2 of the replicates.

```python
sig_per_exp = []
for exp in ace["expName"]:
    df = condition_top_genes(exp, n=200, sign="neg")
    sig_per_exp.append(set(df["locusId"]))

from collections import Counter
counts = Counter()
for s in sig_per_exp:
    counts.update(s)
robust = [loc for loc, c in counts.items() if c >= 2]
print(f"{len(robust)} genes significant in ≥2 acetate replicates")
```

---

## 4. Compare two conditions (or substrate vs control)

> "Which genes matter on acetate but not glucose?"

```python
# pick representative experiments (one each — average reps if you want)
ace = experiments(contains="acetate").iloc[0]["expName"]   # set1IT084
glu = experiments(contains="glucose").iloc[0]["expName"]   # set1IT078

cmp = condition_compare(ace, glu)

# genes specifically required on acetate (more depleted on acetate, fit_a < fit_b, large negative delta)
acetate_specific = cmp.loc[(cmp["fit_a"] < -1) & (cmp["fit_b"] > -0.5)].head(20)
print(acetate_specific.to_string(index=False))
```

---

## 5. Pathway sanity check

> "Does the leucine biosynthesis pathway show up as expected in MOPS-no-Leu?"

```python
leu_pathway = ["PP_3540", "PP_3541", "PP_4807", "PP_4808", "PP_4809"]   # leuABCD etc
mops_minus_leu = experiments(contains="leu")  # filters expDesc; refine as needed
grid = substrate_grid(genes=leu_pathway, substrates=mops_minus_leu["expDesc"].tolist())
plot_substrate_heatmap(grid, out="/tmp/leu_pathway.png")
```

If all genes show fit ≪ 0 in the leu-dropout experiment but ≈ 0 elsewhere → pathway is functional and rate-limiting in this context. ✓

---

## 6. Cofitness-driven hypothesis generation

> "Find genes that move together with PP_0154 — likely co-functional."

```python
neighbors = cofitness("PP_0154", n=15)
print(neighbors.to_string(index=False))

# Expand: get the gene's annotation for each
for nb in neighbors["hitId"]:
    g = gene(nb)
    print(f"  {nb}  {g.get('reannotation', g.get('desc', '?'))[:80]}")
```

`conserved=TRUE` rows are highest-confidence functional partners.

---

## 7. QC a suspect experiment

```python
ok, details = qc_pass("set10IT006")
print(f"PASS: {ok}")
for metric, (val, thresh, passed) in details.items():
    mark = "✓" if passed else "✗"
    print(f"  {mark}  {metric}: {val} (threshold {thresh})")
```

If an experiment fails QC, treat its phenotypes with skepticism — the strongest hits often survive but borderline calls (|fit|=1, |t|=4) shouldn't be reported.

---

## 8. Cross-reference with UniProt / GO

`gene()` already pulls UniProt fields when present:

```python
g = gene("PP_0154")
print(g.get("EC number"))
print(g.get("Pathway"))
print(g.get("Gene Ontology (biological process)"))
```

Use these to call out a phenotype's biological context (e.g. "GO:0006083 acetate metabolism — fits expectation").

---

## Output style

Every analysis should end with:

```
PP_NNNN — <one-line annotation>
  · <strongest specific phenotype, with expDesc + fit + |t|>
  · <top cofit neighbor and what it implies>
  · <link to fitness browser if user might want to explore: https://fit.genomics.lbl.gov/cgi-bin/myFitShow.cgi?orgId=Putida&gene=PP_NNNN>
```

Three lines, no fluff.
