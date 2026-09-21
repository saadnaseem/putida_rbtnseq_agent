# Data provenance and licence

The `data/` directory is redistributed here under the terms below. **The MIT licence in
`LICENSE` covers the code only — it does not apply to these files.**

If you use this data, cite the sources, not this repository. Full references are in
[README.md](README.md#citing).

---

## 1. Fitness Browser files — CC BY 4.0

`fit_organism_Putida.tsv`, `t_organism_Putida.tsv`, `exp_organism_Putida.txt`,
`organism_Putida_genes.tab`, `reanno_Putida.tsv`, `cofit_organism_Putida.txt`,
`specific_phenotypes_Putida.txt`, `organism_Putida.faa`, `organism_Putida.fna`

- **Source:** LBNL Fitness Browser, <https://fit.genomics.lbl.gov/cgi-bin/org.cgi?orgId=Putida>
- **Authors:** Morgan N. Price, Adam M. Deutschbauer, Adam P. Arkin (Lawrence Berkeley
  National Laboratory)
- **Licence:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/), as declared on
  the Fitness Browser's figshare releases (e.g.
  [10.6084/m9.figshare.5134840](https://figshare.com/articles/dataset/5134840))
- **Cite:** Price MN, Wetmore KM, Waters RJ, et al. *Mutant phenotypes for thousands of
  bacterial genes of unknown function.* Nature 557, 503–509 (2018).
  <https://doi.org/10.1038/s41586-018-0124-0>
- **Retrieved:** 21 September 2026, via `scripts/fetch_data.sh`

## 2. fModule supplement — CC BY 4.0

`fModule_Metadata.xlsx`

- **Source:** <https://github.com/beckham-lab/fModule> (the repository named in the
  paper's data availability statement), byte-content verified identical across all three
  sheets to the published supplement
- **Licence:** the article is open access under
  [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/); the hosting repository
  declares GPL-3.0
- **Cite:** Borchert AJ, Bleem AC, Lim HG, Rychel K, Dooley KD, Kellermyer ZA, Hodges TL,
  Palsson BO, Beckham GT. *Machine learning analysis of RB-TnSeq fitness data predicts
  functional gene modules in Pseudomonas putida KT2440.* mSystems 9(3) (2024).
  <https://doi.org/10.1128/msystems.00942-23>
- **Retrieved:** 21 September 2026

## 3. UniProt proteome — CC BY 4.0

`uniprot_putida.csv`

- **Source:** UniProt proteome [UP000000556](https://www.uniprot.org/proteomes/UP000000556)
  (*Pseudomonas putida* KT2440)
- **Licence:** [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)
- **Cite:** The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase.*
  Nucleic Acids Research.

---

## Not redistributed here

The Borchert et al. article PDF is not included. It is open access — read it at the DOI
above.

## A note on personal data

`exp_organism_Putida.txt` carries `person` and `dateStarted` columns naming the
researchers who ran each experiment. These are part of the public Fitness Browser
release and are redistributed unmodified. No code in this repository reads them.

## Refreshing

`scripts/fetch_data.sh` re-downloads the Fitness Browser files. The numeric matrices have
been stable across releases (verified bit-identical between May and September 2026); the
curation layers grow. See [docs/DATA_CAVEATS.md](docs/DATA_CAVEATS.md) §4.
