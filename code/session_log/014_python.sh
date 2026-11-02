# session cell 14 | 2026-10-15 19:48 | env: python | language: bash
# files written: []

curl -s -m 20 "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=nuccore&term=Psychromicrobium%5BOrganism%5D+AND+16S+ribosomal+RNA%5BTitle%5D&retmax=20" | head -c 800