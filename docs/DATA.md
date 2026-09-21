# The dataset: what ships, and how to refresh it

**The data is in this repository.** `data/` holds all 11 source files (71 MB), committed.
After `git clone` you need one command and no network:

```bash
python scripts/build_cache.py    # ~2 min
```

That converts the raw files into the parquet cache (`data/cache/`, 45 MB) that the API
reads. The cache is gitignored because `build_cache.py` regenerates it byte-identically.

Licence and per-file provenance: [DATA_LICENSE.md](../DATA_LICENSE.md). All three sources
are CC BY 4.0.

---

## What's in `data/`

### From the LBNL Fitness Browser

Retrieved 21 Sept 2026 from <https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=Putida>.

| File | Source endpoint (`.../cgi-bin/`) | Contents |
|---|---|---|
| `fit_organism_Putida.tsv` | `createFitData.cgi?orgId=Putida` | Fitness matrix, 4,778 genes × 314 experiments |
| `t_organism_Putida.tsv` | `createFitData.cgi?orgId=Putida&t=1` | Matching t-like statistics |
| `exp_organism_Putida.txt` | `createExpData.cgi?orgId=Putida` | Experiment metadata **+ the QC fields** |
| `organism_Putida_genes.tab` | `orgGenes.cgi?orgId=Putida` | Gene table, 5,661 genes |
| `reanno_Putida.tsv` | `downloadReanno.cgi?orgId=Putida` | 41 curated re-annotations |
| `cofit_organism_Putida.txt` | `createCofitData.cgi?orgId=Putida` | Precomputed top cofitness per gene |
| `specific_phenotypes_Putida.txt` | `spec.cgi?orgId=Putida&download=1` | The browser's own significant calls |
| `organism_Putida.fna` | `orgSeqs.cgi?orgId=Putida&type=nt` | Genome assembly (FASTA) |
| `organism_Putida.faa` | `orgSeqs.cgi?orgId=Putida` | Protein sequences (FASTA) |

The FASTA files are not used by the code; they are included so `data/` is a complete
mirror of what `references/files.md` documents.

> The Fitness Browser now lists this organism under its reclassified name,
> **_Aquipseudomonas alloputida_ KT2440**. The `orgId=Putida` handle is unchanged and it
> is the same strain.

### From the Borchert et al. 2024 supplement

`fModule_Metadata.xlsx`, taken from <https://github.com/beckham-lab/fModule> — the
repository named in the paper's own data availability statement. Its content was verified
identical, sheet by sheet, to the published supplement.

> Borchert AJ, Bleem AC, Lim HG, Rychel K, Dooley KD, Kellermyer ZA, Hodges TL,
> Palsson BO, Beckham GT. **Machine learning analysis of RB-TnSeq fitness data predicts
> functional gene modules in *Pseudomonas putida* KT2440.** *mSystems* 9(3) (2024).
> <https://doi.org/10.1128/msystems.00942-23>

Three sheets, named exactly:

| Sheet | Shape | Contents |
|---|---|---|
| `fitness_measurements` | 4,732 × 337 | Fitness matrix, **332** experiments |
| `T-like_statistics` | 4,732 × 337 | Matching t statistics |
| `metadata` | 332 × 29 | Experiment metadata (no QC fields) |

This is the **canonical** matrix for analysis — a superset of the public snapshot
(332 vs 314 experiments). Read [DATA_CAVEATS.md](DATA_CAVEATS.md) before interpreting its
`t` column.

### From UniProt

`uniprot_putida.csv` — proteome [UP000000556](https://www.uniprot.org/proteomes/UP000000556),
used to attach EC numbers, GO terms, pathways and subcellular location. The join key is
`Gene Names (ordered locus)` → `locusId`; `build_cache.py` renames it.

---

## Refreshing from upstream

```bash
./scripts/fetch_data.sh          # re-downloads the 7 Fitness Browser tables into data/
python scripts/build_cache.py --rebuild
pytest -q
```

The script skips files that are already present, so delete the ones you want replaced
first. It does not fetch the FASTA, the xlsx or the UniProt export — those change rarely;
the table above and the sections above give their URLs.

**What changes between releases.** Verified by rebuilding from scratch in September 2026
against a May 2026 copy: all 3,071,316 fitness and t values came back **bit-identical**,
as did the gene table, metadata and UniProt export. Only curation layers moved —
`specific_phenotypes` +37 calls (none removed), two `conserved` flags, and the organism
rename. See [DATA_CAVEATS.md](DATA_CAVEATS.md) §4.

`data/cache/MANIFEST.json` records row/column counts and a checksum per cached file, so
you can pin exactly what an analysis ran against.

---

## Using a dataset stored elsewhere

The code resolves its data directory as `$PUTIDA_RBTNSEQ_DATA` if set, otherwise
`<repo>/data`. To point at a copy outside the repo:

```bash
export PUTIDA_RBTNSEQ_DATA=~/data/rbtnseq
```

---

## Verify

```bash
pytest -q
```

25 tests, covering matrix dimensions, join integrity, the documented data caveats, and
that the significance rule reproduces the Fitness Browser's own specific-phenotype calls.
If the cache is missing, the data-dependent tests skip rather than fail.

---

## Note on personal data

`exp_organism_Putida.txt` carries `person` and `dateStarted` columns naming the
researchers who ran each experiment. These are part of the public Fitness Browser release
and are redistributed unmodified. No code here reads them.
