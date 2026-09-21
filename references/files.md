# Data file inventory

Raw data lives at `$PUTIDA_RBTNSEQ_DATA/`. Parquet cache (auto-built by `scripts/build_cache.py`) lives at `$PUTIDA_RBTNSEQ_DATA/cache/`.

**You almost never need to read raw files.** Use `scripts/putida_rbtnseq.py` — it reads the parquet cache.

## Parquet cache (load via helper module)

| Cache file | Source | Used by |
|---|---|---|
| `genes.parquet` | `organism_Putida_genes.tab` | `gene()` |
| `reanno.parquet` | `reanno_Putida.tsv` | `gene()` (overlay) |
| `uniprot.parquet` | `Pseudomonas_putida_protein_list_160488.csv` | `gene()` (overlay) |
| `metadata_paper.parquet` | `fModule_Metadata.xlsx → metadata` | `experiments()`, `phenotypes()` (paper, 332 expts) |
| `metadata.parquet` | `exp_organism_Putida.txt` | `qc_pass()`, fallback (public, 314 expts) |
| `fit_paper.parquet` / `t_paper.parquet` | xlsx sheets | wide matrices (genes × expts) |
| `fit_paper_long.parquet` | derived | tidy long form (locusId, expName, fit, t) — what most queries hit |
| `fit.parquet` / `t.parquet` / `fit_long.parquet` | public TSVs | same shape, public-snapshot version |
| `cofit.parquet` | `cofit_organism_Putida.txt` | `cofitness()` |
| `specific_phenotypes.parquet` | `specific_phenotypes_Putida.txt` | optional pre-filtered list of significant calls |

## Raw files (if you really need them)

### Genome / annotation
- `organism_Putida.fna` — genome assembly (1 scaffold, `AE015451`)
- `organism_Putida.faa` — 5,563 protein FASTAs, headers `>Putida:PP_NNNN  PP_NNNN  desc`
- `organism_Putida_genes.tab` — 5,661 gene table: `locusId, sysName, scaffoldId, begin, end, strand, type, desc`
- `reanno_Putida.tsv` — 41 curated re-annotations: `locusId, gene, reannotation, comment, original_desc, aaseq`
- `Pseudomonas_putida_protein_list_160488.csv` — 5,528 UniProt rows: EC, GO, pathway, function, subcellular, transmembrane, etc. Join on `Gene Names (ordered locus)` ↔ `locusId`

### Experimental metadata
- `exp_organism_Putida.txt` — 314 experiments. **QC fields live here, not in the paper xlsx**: `gMed`, `mad12`, `cor12`, `gccor`, `adjcor` (plus `mad12c`, `opcor`, which are *not* QC gates — see biology.md), `nMapped`, `nUsed`. Other cols: `expName`, `expDesc`, `expGroup`, `media`, `temperature`, `pH`, `vessel`, `aerobic`, `liquid`, `condition_1..4`, `units_1..4`, `concentration_1..4`, `growthPlate`, `growthWells`
- `fModule_Metadata.xlsx → metadata` — 332 experiments (paper version). Adds `set`, `classifier`, `total_rep`, `rep`, `Inoculum media type`. **No QC fields.**

### Fitness matrices
- `fit_organism_Putida (1).tsv` — 4,778 genes × 314 experiments. First 5 cols `orgId, locusId, sysName, geneName, desc`; remaining cols are `<expName> <expDesc>` (space-separated)
- `t_organism_Putida.tsv` — same shape, t-statistics
- `fModule_Metadata.xlsx → fitness_measurements` / `T-like_statistics` — paper version, 4,732 × 332

### Derived
- `cofit_organism_Putida.txt` — 95,560 rows. Per-gene top cofit neighbors (Pearson corr of fitness profiles). Cols: `orgId, locusId, sysName, name, desc, hitId, hitSysName, hitName, hitDesc, rank, cofit, conserved`. `conserved=TRUE` = same cofit relationship in orthologs (high-confidence).
- `specific_phenotypes_Putida.txt` — 2,437 pre-filtered (|fit|>1, |t|>4) gene-condition phenotypes: `expGroup, expName, condition_1, locusId, sysName, gene, desc, fit, t, conserved`

### Background
- `f-gene-modules.pdf` — Borchert et al. 2024, *mSystems*: ICA on this fitness matrix → 84 fModules. The xlsx is the input matrix.

## Join keys

```
genes / reanno / uniprot / fit_long / cofit / specific_phenotypes
                ↑
            locusId  (e.g. PP_0154)

metadata / fit_long / specific_phenotypes
                ↑
            expName  (e.g. set10IT006)
```

## Rebuild cache

```bash
python \
    $SKILL/scripts/build_cache.py --rebuild
```

Takes ~45 s (xlsx parse is the bottleneck). Idempotent without `--rebuild`.
