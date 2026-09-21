# RB-TnSeq biology primer & interpretation

## What's measured

A pooled library of ~10⁵ *P. putida* mutants (each with one transposon insertion + a 20-bp barcode) is grown competitively in a defined condition. Barcodes are PCR-amplified and sequenced at:

- **t₀** — inoculum / time-zero sample
- **t_end** — end of growth (typically 4–8 generations)

For each gene, fitness is computed across all insertion strains in that gene:

```
fit(gene, condition) ≈ log2( median end-count / median t0-count )    (per-strain, then aggregated)
```

The values are normalized so:
- the chromosomal trend across position is removed (`gccor`)
- the median essential-gene set sits at fitness ≈ 0
- a typical non-essential gene also sits near 0

So a gene with `fit ≈ 0` is just neutral. A gene with `fit < -2` is strongly required. A gene with `fit > +2` is strongly inhibitory in that condition (loss helps).

## What "t" means

`t` is a t-like statistic that asks: is this gene's fitness reliably different from zero, given the within-gene noise (variance across the multiple insertion strains in that gene) and read depth?

- |t| < 2  →  noise
- |t| ~ 4  →  reasonably confident
- |t| > 4  →  Arkin lab's "specific phenotype" cutoff

**Unsigned `t` in two experiment sets**: in the paper xlsx
(`fModule_Metadata.xlsx → T-like_statistics`), **`set100` and `set101` store `t` as an
unsigned magnitude**. Verified: zero negative values across all 369,096 measurements in
those two sets, while every other set is ~52% negative. This affects **78 of the 332
paper experiments (23.5%)** — none of which appear in the public snapshot, which is
100% sign-consistent.

So `set100IT039` reads `fit = −1.603, t = +7.331`: the gene is depleted (negative `fit`)
with high confidence (`|t|` = 7.3). The positive sign on `t` carries no meaning.

The helper functions threshold on **`|t|`**, which is correct for both cases.
**Do not interpret signed `t` directly** — use signed `fit` for direction and `|t|` for
confidence. Call `has_signed_t(expName)` to test an experiment, or see
`SETS_WITHOUT_SIGNED_T`. A regression test pins this invariant.

## Standard significance cutoff

A gene is a "specific phenotype" in a condition when:

```
|fit| > 1  AND  |t| > 4
```

This is what the Fitness Browser uses to populate its "specific phenotypes" tables and what `specific_phenotypes_Putida.txt` is pre-filtered on.

For looser exploratory work, |fit|>0.5 & |t|>3 is sometimes used; for stringent work, |fit|>2 & |t|>5.

## Interpretation cheat sheet

| Pattern | Likely interpretation |
|---|---|
| fit ≪ 0, |t| > 4 | Gene is **important for growth** in this condition. Knockout is depleted. Could be biosynthesis (made by gene → media lacks it), uptake, central metabolism, stress resistance |
| fit ≫ 0, |t| > 4 | Gene **inhibits growth** in this condition. Knockout is enriched. Often regulators, repressors, or genes whose products are toxic in that condition |
| fit ≈ 0 across all conditions | Gene is dispensable / redundant / not expressed |
| fit ≈ 0 in test, but very low in t₀ control | Likely **essential** gene — too few insertions to score in any condition. Look at `nMapped` for that gene |
| fit < 0 in many conditions | Generally important / housekeeping |
| fit < 0 in *only* one substrate condition | Specific catabolic / uptake gene for that substrate |

## Experiment QC thresholds (from `exp_organism_Putida.txt`)

A trustworthy experiment passes:

| Metric | Threshold | Meaning |
|---|---|---|
| `gMed` | ≥ 50 | Median reads per gene — adequate sequencing depth |
| `mad12` | ≤ 0.5 | Median absolute difference in fitness between the two halves of each gene |
| `cor12` | ≥ 0.1 | Correlation in fitness between the two halves of each gene |
| `gccor` | \|gccor\| ≤ 0.2 | Correlation between gene GC content and fitness — bias check |
| `adjcor` | \|adjcor\| ≤ 0.25 | Correlation between adjacent genes — bias check |

These are the criteria from Wetmore et al. 2015 (mBio), as stated in the Fitness Browser help.

> **Use `mad12`, not `mad12c`.** Despite the name, `mad12c` is a correlation-like
> quantity centred near 0.96 (range 0.91–1.10 across all 314 experiments), not a
> median absolute difference. Testing it against `≤ 0.5` rejects 100% of experiments.
> `mad12` is the real metric (centred 0.20, max 0.46). There is likewise no published
> `opcor − adjcor ≥ 0.1` criterion; `opcor`/`adjcor` are reported for information.

**The public snapshot is already QC-filtered.** `gMed` bottoms out at exactly 50.0 and
`cor12` at exactly 0.101 — the Fitness Browser applies these cutoffs before release, so
all 314 experiments pass by construction. `qc_pass()` is a guard for re-derived or newly
added data, not a filter that will reject anything in the shipped snapshot.

Use `qc_pass(expName)` in the helper module.

The **paper xlsx does not ship QC fields**, so QC always uses the public-snapshot metadata. Most experiments overlap between the two, but a few paper-only experiments can't be QC'd.

## Cofitness

`cofit(gene_a, gene_b)` = Pearson correlation of their fitness vectors across all 314 experiments.

- cofit > 0.7 → very likely same pathway / same complex / regulon
- cofit 0.4–0.7 → functionally related (same general process)
- cofit < 0.3 → uninformative (any random pair sits here)
- `conserved=TRUE` → same relationship seen in orthologs in other organisms — much higher confidence

The Fitness Browser ships pre-computed top-N cofit per gene in `cofit_organism_Putida.txt`. Use `cofitness(locus, n=20)` to read it. (Computed on the 314-experiment public snapshot; not recomputed for the paper's 332-experiment matrix.)

## fModules (Borchert 2024)

Independent Component Analysis (ICA) was run on the fitness matrix. Each independent component → an "fModule" = a group of genes whose fitness moves together. 84 modules total. fModules generalize cofitness from pairwise to multi-gene groups.

We don't ship the fModule assignments as a table — only the input matrices (`fit_paper.parquet`, `t_paper.parquet`). To use modules, re-derive with sklearn's FastICA on `fit_paper.parquet`, or use `cofitness()` as a pairwise proxy. See `fmodules.md`.

## Common pitfalls

1. **Polar effects**: insertions in upstream genes can disrupt downstream genes in the same operon, inflating apparent fitness of the downstream gene. Always check operon context (`organism_Putida_genes.tab` ordering) before claiming a single-gene phenotype.
2. **Essential genes**: have few/no insertions, so they get **NaN** or near-zero scores everywhere — they look "uninformative" even when critical. Look at the gene's `nUsed` in the underlying mapping if a key gene is unexpectedly silent.
3. **Replicates**: many conditions have 2–4 replicates (different `expName`s, same `expDesc`). Average them when reporting; don't double-count.
4. **t₀ pairing**: each experiment is normalized to its own time-zero set (`timeZeroSet`). Cross-set comparisons are valid because of the chromosomal-trend normalization, but very small fitness differences (<0.5) between sets can be batch noise.
