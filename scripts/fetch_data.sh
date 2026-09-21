#!/usr/bin/env bash
# Download the public P. putida KT2440 RB-TnSeq snapshot from the LBNL Fitness Browser.
#
# Usage:
#   ./scripts/fetch_data.sh [target_dir]      # default: $PUTIDA_RBTNSEQ_DATA, else ./data
#
# Downloads the 7 Fitness Browser files (~35 MB). Two further files are NOT
# downloadable this way and must be obtained manually — see docs/DATA.md:
#   * fModule_Metadata.xlsx                     (Borchert et al. 2024 supplement)
#   * Pseudomonas_putida_protein_list_*.csv     (UniProt proteome export)

set -euo pipefail

DEST="${1:-${PUTIDA_RBTNSEQ_DATA:-./data}}"
BASE="https://fit.genomics.lbl.gov/cgi-bin"
ORG="Putida"

mkdir -p "$DEST"
cd "$DEST"
echo "Downloading into: $(pwd)"

fetch () {  # fetch <url-path> <output-name>
  if [[ -s "$2" ]]; then
    echo "  [skip] $2 (already present)"
  else
    echo "  [get ] $2"
    curl -fsSL --retry 3 "$BASE/$1" -o "$2"
  fi
}

fetch "createFitData.cgi?orgId=$ORG"          "fit_organism_${ORG}.tsv"
fetch "createFitData.cgi?orgId=$ORG&t=1"      "t_organism_${ORG}.tsv"
fetch "createExpData.cgi?orgId=$ORG"          "exp_organism_${ORG}.txt"
fetch "orgGenes.cgi?orgId=$ORG"               "organism_${ORG}_genes.tab"
fetch "downloadReanno.cgi?orgId=$ORG"         "reanno_${ORG}.tsv"
fetch "createCofitData.cgi?orgId=$ORG"        "cofit_organism_${ORG}.txt"
fetch "spec.cgi?orgId=$ORG&download=1"        "specific_phenotypes_${ORG}.txt"

echo
echo "Done. Still required (manual — see docs/DATA.md):"
[[ -f fModule_Metadata.xlsx ]] \
  && echo "  [ok] fModule_Metadata.xlsx" \
  || echo "  [ MISSING ] fModule_Metadata.xlsx  (Borchert 2024 supplementary data)"
ls Pseudomonas_putida_protein_list_*.csv uniprot_putida.csv >/dev/null 2>&1 \
  && echo "  [ok] UniProt proteome CSV" \
  || echo "  [ MISSING ] UniProt proteome CSV  (uniprot.org export)"
echo
echo "Then build the cache:  python scripts/build_cache.py"
