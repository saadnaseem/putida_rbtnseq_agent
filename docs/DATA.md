# Obtaining the dataset

**No data is redistributed in this repository.** The sources have their own terms and
one of them (a journal supplement) is not ours to mirror. This document is the complete
protocol for reconstructing the exact dataset the code expects.

Total download: **~120 MB**. Build time: **~2 minutes**.

---

## 1. Choose where the data lives

```bash
export PUTIDA_RBTNSEQ_DATA=~/data/rbtnseq     # anywhere you like
```

If you skip this, the code defaults to `<repo>/data`, which is gitignored.

---

## 2. Fitness Browser files (automated, ~35 MB)

```bash
./scripts/fetch_data.sh
```

This pulls seven files from the LBNL Fitness Browser
(<https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=Putida>):

| Output file | Source endpoint (`.../cgi-bin/`) | Contents |
|---|---|---|
| `fit_organism_Putida.tsv` | `createFitData.cgi?orgId=Putida` | Fitness matrix, 4,778 genes × 314 experiments |
| `t_organism_Putida.tsv` | `createFitData.cgi?orgId=Putida&t=1` | Matching t-like statistics |
| `exp_organism_Putida.txt` | `createExpData.cgi?orgId=Putida` | Experiment metadata **+ the QC fields** |
| `organism_Putida_genes.tab` | `orgGenes.cgi?orgId=Putida` | Gene table, 5,661 genes |
| `reanno_Putida.tsv` | `downloadReanno.cgi?orgId=Putida` | 41 curated re-annotations |
| `cofit_organism_Putida.txt` | `createCofitData.cgi?orgId=Putida` | Precomputed top cofitness per gene |
| `specific_phenotypes_Putida.txt` | `spec.cgi?orgId=Putida&download=1` | The browser's own significant calls |

> The Fitness Browser now lists this organism under its reclassified name,
> **_Aquipseudomonas alloputida_ KT2440**. The `orgId=Putida` handle is unchanged
> and it is the same strain.

---

## 3. Borchert et al. 2024 supplement (manual, ~23 MB)

`fModule_Metadata.xlsx` is supplementary data from:

> Borchert AJ, Bleem AC, Lim HG, et al. **Machine learning and systems biology
> approaches reveal fitness modules in *Pseudomonas putida* KT2440.**
> *mSystems* (2024). <https://doi.org/10.1128/msystems.00934-24>

Download the supplementary spreadsheet from the article page and save it as
`$PUTIDA_RBTNSEQ_DATA/fModule_Metadata.xlsx`.

It must contain three sheets, named exactly:

| Sheet | Shape | Contents |
|---|---|---|
| `fitness_measurements` | 4,732 × 337 | Fitness matrix, **332** experiments |
| `T-like_statistics` | 4,732 × 337 | Matching t statistics |
| `metadata` | 332 × 29 | Experiment metadata (no QC fields) |

This is the **canonical** matrix for analysis — it is a superset of the public
snapshot (332 vs 314 experiments). See [DATA_CAVEATS.md](DATA_CAVEATS.md) before
interpreting its `t` column.

---

## 4. UniProt proteome (manual, ~6 MB)

Used to attach EC numbers, GO terms, pathways and subcellular location to genes.

1. Go to <https://www.uniprot.org/proteomes/UP000000556> (*P. putida* KT2440).
2. View all proteins, then **Download → Format: CSV → Uncompressed**.
3. Include at least these columns: `Entry`, `Protein names`, `Gene Names (ordered locus)`,
   `EC number`, `Pathway`, `Function [CC]`, the three `Gene Ontology` columns, and
   `Subcellular location [CC]`.
4. Save as `$PUTIDA_RBTNSEQ_DATA/uniprot_putida.csv`
   (a file named `Pseudomonas_putida_protein_list_*.csv` is also recognised).

The join key is `Gene Names (ordered locus)` → `locusId`; `build_cache.py` renames it.

---

## 5. Build the cache

```bash
python scripts/build_cache.py
```

Expected output (~2 min, dominated by the xlsx read):

```
genes.parquet                          0.22 MB
reanno.parquet                         0.03 MB
uniprot.parquet                        2.83 MB
metadata.parquet                       0.07 MB
fit.parquet                            4.11 MB
t.parquet                              6.56 MB
fit_long.parquet                      10.38 MB
fit_paper.parquet                      4.19 MB
t_paper.parquet                        6.50 MB
fit_paper_long.parquet                10.49 MB
metadata_paper.parquet                 0.03 MB
cofit.parquet                          1.50 MB
specific_phenotypes.parquet            0.08 MB
```

`cache/MANIFEST.json` records row/column counts and a checksum per file so a rebuild
can be compared against a known-good one.

---

## 6. Verify

```bash
pytest -q
```

24 tests. They check matrix dimensions, join integrity, the documented data caveats, and
that the significance rule reproduces the Fitness Browser's own specific-phenotype calls.
If the cache is absent the data-dependent tests skip rather than fail.

---

## Note on personal data

`exp_organism_Putida.txt` includes a `person` column naming the researcher who ran each
experiment, and `dateStarted`. These are part of the public Fitness Browser release. This
repository neither redistributes nor uses those columns, but be aware they are present in
your local copy if you plan to share it.
