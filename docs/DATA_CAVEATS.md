# Data caveats

Two properties of the source data will silently produce wrong answers if you are not
aware of them. Both were found by auditing this dataset; both are pinned by regression
tests in `tests/test_putida_rbtnseq.py`.

---

## 1. `set100` and `set101` store `t` as an unsigned magnitude

### What you see

```
PP_0154  set100IT039  "D-Glucose with Acetic acid (C)"   fit = -1.603   t = +7.331
```

A negative fitness with a positive t-statistic is not possible from a single
measurement — `t` is the fitness divided by its standard error, so the two must share
a sign.

### What is actually happening

In `fModule_Metadata.xlsx → T-like_statistics`, two experiment sets report `t` without
its sign:

| Set | measurements | negative `t` values | min `t` |
|---|---|---|---|
| `set100` | 184,548 | **0** | 0.000 |
| `set101` | 184,548 | **0** | 0.000 |
| every other set | — | ~52% | ≈ −20 |

Zero negative values out of 369,096 is not a coincidence; those two sets carry `|t|`.
The magnitudes themselves are sound — within `set100`, `|fit|` and `|t|` correlate at
0.73, essentially the same as the 0.78 seen in a well-behaved set.

**Scope:** 78 of the 332 paper experiments (23.5%). None of them appear in the public
Fitness Browser snapshot, which is 100% sign-consistent across all 27,013 of its
significant calls.

### What to do

- Use signed **`fit`** for direction (negative = mutant depleted = gene required).
- Use **`|t|`** for confidence only.
- Never branch on `sign(t)`.

All helper functions already threshold on `|t|`, so `phenotypes()`,
`condition_top_genes()` and `condition_compare()` are correct as written. To test an
experiment directly:

```python
>>> has_signed_t("set1IT084")
True
>>> has_signed_t("set100IT039")
False
```

Guarded by `test_set100_and_set101_store_unsigned_t` and
`test_public_snapshot_is_sign_consistent`.

---

## 2. QC uses `mad12`, not `mad12c`

### The trap

The published criterion is *"the median absolute difference in fitness between the two
halves of a gene is ≤ 0.5"*. Two columns look like candidates. Only one is right:

| Column | min | median | max | usable as a `≤ 0.5` gate? |
|---|---|---|---|---|
| `mad12` | 0.144 | 0.196 | 0.458 | **yes** |
| `mad12c` | 0.914 | 0.956 | 1.100 | no — rejects 100% of experiments |

Despite the name, `mad12c` is a correlation-like quantity, not a median absolute
difference. An earlier version of this code tested `mad12c ≤ 0.5` and therefore failed
**314 of 314** experiments.

There is also no published `opcor − adjcor ≥ 0.1` criterion. That rule rejects a further
27.7% of experiments for no documented reason. `opcor` and `adjcor` are reported for
information; only `|adjcor| ≤ 0.25` is an actual gate.

### The criteria actually used

From Wetmore et al. 2015 (*mBio*), as stated in the Fitness Browser help:

| Metric | Gate | Meaning |
|---|---|---|
| `gMed` | ≥ 50 | median reads per gene |
| `mad12` | ≤ 0.5 | median absolute fitness difference between gene halves |
| `cor12` | ≥ 0.1 | correlation between gene halves |
| `gccor` | \|·\| ≤ 0.2 | fitness-vs-GC-content bias |
| `adjcor` | \|·\| ≤ 0.25 | adjacent-gene correlation |

All 314 public experiments pass.

### Why they all pass

The Fitness Browser applies these cutoffs **before release**: `gMed` bottoms out at
exactly 50.0 and `cor12` at exactly 0.101. So `qc_pass()` will not reject anything in
the shipped snapshot. It is a guard for re-derived or newly added experiments, not a
filter to run over the distribution.

The 78 paper-only experiments cannot be QC'd at all — the paper xlsx ships no QC fields.
`qc_pass()` raises a `KeyError` that says so.

Guarded by `test_qc_metrics_use_documented_columns`,
`test_every_public_experiment_passes_qc` and `test_qc_on_paper_only_experiment_explains_itself`.

---

## 3. Smaller things worth knowing

- **929 of 5,661 annotated genes have no fitness data** — essential genes, or genes with
  too few transposon insertions to score. `gene()` returns their annotation with
  `has_fitness_data=False` rather than raising. Absence of a phenotype here is a
  measurement gap, not evidence of dispensability.
- **The two matrices differ in both dimensions**: paper 4,732 × 332, public 4,778 × 314.
  Do not join them positionally; join on `locusId` / `expName`.
- **Cofitness is computed on the 314-experiment basis** even when you are analysing the
  332-experiment paper matrix. The paper ships no recomputed cofitness table.
- **Replicates share an `expDesc`** (e.g. `set1IT084` and `set1IT085` are both
  "Acetate Carbon Source (20mM)"). Average them before reporting; don't double-count.
- **Polar effects**: an insertion can disrupt downstream genes in the same operon. Check
  operon context before claiming a single-gene phenotype.
