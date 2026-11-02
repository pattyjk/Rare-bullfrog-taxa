# session cell 162 | 2026-10-18 12:52 | env: r | language: bash
# files written: []

cd amphib_R; grep -n "v4_window\|hamming_to\|norm_bd\|MIN_READS\|<- function\|FWD\|primer" R/00_functions.R | head -30; sed -n '/v4_window <- /,/^}/p;/norm_bd <- /,/^}/p' R/00_functions.R; grep -n "3mm\|mism\|containment\|grepl\|stringdist\|inhib_score" R/04_amphibac_match.R | head -20; wc -l < data/derived/sample_totals.csv