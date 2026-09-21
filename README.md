# putida_rbtnseq_agent

Query and analyse **RB-TnSeq gene-fitness data for _Pseudomonas putida_ KT2440** from the
LBNL Fitness Browser and the Borchert et al. 2024 fModules paper — as a Python API and as
a [Claude Code](https://claude.com/claude-code) skill.

Ask *"what does PP_0154 do?"* and get an answer grounded in 332 competitive-fitness
experiments instead of a guess from the gene name.

```
PP_0154 — propionyl-CoA:succinate CoA transferase
  · Specifically required for acetate growth (set1IT084: fit=-1.52, |t|=9.9)
  · Top cofit: PP_4487 acetyl-CoA synthetase (cofit=0.57) → acetate-uptake module
```

---

## What this is

**RB-TnSeq** (Random-Barcoded Transposon Sequencing) measures, for every gene in the
genome at once, how much a knockout mutant is depleted or enriched when the population is
grown competitively under a defined condition. A pooled library of ~10⁵ barcoded
insertion mutants is grown, barcodes are counted by sequencing at the start and end, and
each gene gets:

- **`fit`** — log₂(end / start), normalised so a typical gene sits at 0.
  Negative = mutant depleted = **gene was needed**. Positive = its loss helped.
- **`t`** — a t-like statistic: is that fitness reliably different from zero?

This repository gives you fast, correct access to that data, with the dataset's real
quirks handled rather than papered over.

| | |
|---|---|
| Organism | *P. putida* KT2440 (Fitness Browser `orgId=Putida`) |
| Experiments | **332** (paper matrix) / **314** (public snapshot) |
| Genes scored | 4,732 of 5,661 annotated |
| Conditions | carbon sources, nitrogen sources, stresses, bioreactor, supernatants |
| Load time | ~50 ms from the parquet cache |

---

## Quickstart

```bash
git clone https://github.com/saadnaseem/putida_rbtnseq_agent.git
cd putida_rbtnseq_agent
pip install -r requirements.txt

export PUTIDA_RBTNSEQ_DATA=~/data/rbtnseq    # where the dataset will live
./scripts/fetch_data.sh                      # Fitness Browser files (~35 MB)
# two files need a manual download — see docs/DATA.md
python scripts/build_cache.py                # ~2 min
pytest -q                                    # 25 tests
```

Then:

```python
import sys; sys.path.insert(0, "scripts")
from putida_rbtnseq import gene, phenotypes, cofitness

gene("PP_0154")["desc"]          # 'propionyl-CoA:succinate CoA transferase'
phenotypes("PP_0154")            # the 5 conditions where it matters
cofitness("PP_0154", n=5)        # its functional neighbours
```

`python scripts/putida_rbtnseq.py` runs a smoke test and prints exactly this:

```
Smoke test on PP_0154 (propionyl-CoA:succinate CoA transferase)
  desc: propionyl-CoA:succinate CoA transferase
  mean_fit: -0.016, min: -1.603, max: 0.535
  n_phenotypes: 5

Top 5 phenotypes:
    expName                        expDesc      expGroup    fit      t
set100IT039 D-Glucose with Acetic acid (C) carbon source -1.603  7.331
  set1IT084   Acetate Carbon Source (20mM) carbon source -1.521 -9.934
  set1IT085   Acetate Carbon Source (20mM) carbon source -1.479 -8.827
set100IT038 D-Glucose with Acetic acid (C) carbon source -1.414  6.508
set100IT040 D-Glucose with Acetic acid (C) carbon source -1.192  6.816

Top 5 cofit neighbors:
  hitId hitSysName                                hitDesc    cofit  rank  conserved
PP_4487    PP_4487                  acetyl-CoA synthetase 0.568948     1      False
PP_1695    PP_1695  putative Sodium-solute symporter/...  0.432053     2      False
PP_1635    PP_1635         DNA-binding response regulator 0.424029     3      False
PP_1743    PP_1743                       acetate permease 0.366323     4      False
PP_3612    PP_3612       putative TonB-dependent receptor 0.349008     5      False
```

Read that as: knocking out PP_0154 is costly specifically when acetate is the carbon
source, and the genes whose fitness tracks it across all 314 experiments are an acetyl-CoA
synthetase and an acetate permease. That is an acetate-utilisation module, recovered from
data rather than from annotation.

> The `t = +7.331` on `set100IT039` alongside a negative `fit` is **not** a bug —
> `set100`/`set101` store `t` unsigned. See [docs/DATA_CAVEATS.md](docs/DATA_CAVEATS.md).

---

## Data

**No data is redistributed here.** The sources carry their own terms and one is a journal
supplement. [**docs/DATA.md**](docs/DATA.md) is the complete acquisition protocol:
7 files scripted, 2 manual, ~120 MB, ~2 minutes to build.

Data location resolves as `$PUTIDA_RBTNSEQ_DATA`, else `<repo>/data` (gitignored).

---

## API

All analysis goes through `scripts/putida_rbtnseq.py`. **Never read the raw TSV/XLSX** —
the parquet cache loads in under 50 ms.

| Question | Call | Follow up with |
|---|---|---|
| What is this gene, and does it matter anywhere? | `gene("PP_0154")` | `plot_gene_profile(locus, out=...)` |
| Where does it matter? | `phenotypes("PP_0154")` | |
| What moves with it? | `cofitness("PP_0154", n=20)` | |
| Which experiments test substrate X? | `experiments(contains="acetate")` | `experiments(group="carbon source")` |
| What matters in this one condition? | `condition_top_genes("set1IT084", sign="neg")` | `sign="pos"` for genes whose loss helps |
| What changed between two conditions? | `condition_compare(exp_a, exp_b)` | sort by `delta` |
| Genes × substrates map | `substrate_grid(genes, substrates)` | `plot_substrate_heatmap(grid, out=...)` |
| Is this experiment trustworthy? | `qc_pass("set1IT084")` | |
| Is this experiment's `t` signed? | `has_signed_t("set100IT039")` | |

Runnable examples:

```bash
python scripts/examples/gene_profile.py PP_0154
python scripts/examples/condition_phenotypes.py set1IT084
python scripts/examples/substrate_map.py PP_0154,PP_4487,PP_1743 acetate,glucose,citrate
```

### Significance

A gene–condition phenotype is **specific** when `|fit| > 1 AND |t| > 4` — the Arkin lab
convention. Override per call with `fit_thresh=` / `t_thresh=`.

Negative `fit` = mutant depleted = gene important for growth. Use `|t|` for confidence,
never `sign(t)`.

---

## Read before you interpret

[**docs/DATA_CAVEATS.md**](docs/DATA_CAVEATS.md) — two properties of this dataset will
silently produce wrong answers:

1. **`set100` and `set101` store `t` unsigned** (78 of 332 experiments). Verified: zero
   negative values across 369,096 measurements, against ~52% in every other set. Use
   signed `fit` for direction and `|t|` for confidence.
2. **QC uses `mad12`, not `mad12c`.** `mad12c` is a correlation-like value centred at
   0.96; gating it at `≤ 0.5` rejects every experiment. The criteria here are those of
   Wetmore et al. 2015, under which all 314 public experiments pass — as they must, since
   the Fitness Browser pre-filters on them.

Both are pinned by regression tests.

Also in [`references/`](references/): `biology.md` (RB-TnSeq primer, interpretation cheat
sheet, pitfalls), `files.md` (schemas and join keys), `workflows.md` (end-to-end recipes),
`fmodules.md` (Borchert 2024 ICA context).

---

## Validation

```bash
pytest -q     # 25 passed
```

Beyond unit coverage, the suite checks the pipeline against an external ground truth:
applying `|fit| > 1 & |t| > 4` to the public matrix recovers **2,473 of the Fitness
Browser's own 2,474 specific-phenotype calls (99.96%)**. That exercises the whole
load → melt → join path against numbers this repository did not compute. The single miss
is a rounding artifact, not an error — the matrices are rounded to 3 dp while the
phenotype table is not, so a call at `fit = −1.000258` arrives as `−1.000`. The test
asserts that every miss sits within one rounding unit of the thresholds, so a real
pipeline break still fails loudly.

It also pins matrix dimensions, join integrity (no duplicate gene×experiment pairs, no
missing values), the data caveats above, and graceful handling of the 929 genes with no
fitness data.

**Reproducibility was checked end to end**: following `docs/DATA.md` from scratch in
September 2026 against a May 2026 copy reproduced all 3,071,316 fitness and t values
**bit-identically** (max |Δ| = 0). Only the curation layers had grown — `specific_phenotypes`
+37 calls, two `conserved` flags, and the organism rename to *Aquipseudomonas alloputida*.
Details in [docs/DATA_CAVEATS.md](docs/DATA_CAVEATS.md).

Verified on **pandas 2.3.3** (test suite) and **pandas 3.0.3** (full API pass, warnings as
errors). Data-dependent tests skip when the cache is absent, so the suite is safe in CI.

---

## Use as a Claude Code skill

The repository is laid out as a skill — clone it into your skills directory:

```bash
git clone https://github.com/saadnaseem/putida_rbtnseq_agent.git \
  ~/.claude/skills/putida-rbtnseq
```

`SKILL.md` carries the routing rules, a decision tree from question to function call, and
the response format. Claude will load it when you ask about a `PP_` locus, a substrate
phenotype, cofitness, or an fModule. Set `PUTIDA_RBTNSEQ_DATA` in your environment so the
skill can find the cache.

---

## Layout

```
scripts/
  putida_rbtnseq.py      the API — everything goes through here
  build_cache.py         raw sources → parquet cache (+ MANIFEST.json)
  fetch_data.sh          download the Fitness Browser files
  examples/              three runnable end-to-end scripts
references/              biology, file schemas, workflows, fModules
docs/
  DATA.md                acquisition protocol
  DATA_CAVEATS.md        the traps, with evidence
tests/                   24 regression tests
SKILL.md                 Claude Code skill definition
```

---

## Citing

If this is useful, cite the data — not this wrapper:

- **Borchert AJ, Bleem AC, Lim HG, et al.** Machine learning and systems biology approaches
  reveal fitness modules in *Pseudomonas putida* KT2440. *mSystems* (2024).
  <https://doi.org/10.1128/msystems.00934-24>
- **Price MN, Wetmore KM, Waters RJ, et al.** Mutant phenotypes for thousands of bacterial
  genes of unknown function. *Nature* 557, 503–509 (2018).
  <https://doi.org/10.1038/s41586-018-0124-0>
- **Wetmore KM, Price MN, Waters RJ, et al.** Rapid quantification of mutant fitness in
  diverse bacteria by sequencing randomly bar-coded transposons. *mBio* 6:e00306-15 (2015).
  <https://doi.org/10.1128/mBio.00306-15>

Fitness Browser: <https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=Putida>

---

## Licence

Code: [MIT](LICENSE).

The **data is not covered by this licence** and is not distributed here. Fitness Browser
data is produced by the Arkin and Deutschbauer labs at Lawrence Berkeley National
Laboratory; observe the terms on their site and the journal supplements when
redistributing anything derived from them.
