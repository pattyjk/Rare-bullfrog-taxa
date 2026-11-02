# session cell 163 | 2026-10-18 12:52 | env: r | language: bash
# files written: []

cd amphib_R; grep -n "PRIMER_515F <-" R/00_functions.R; sed -n 1,24p R/04_amphibac_match.R | grep -n "exact\|contain" ; Rscript -e 'library(data.table);s<-fread("data/derived/sample_totals.csv");print(names(s));print(nrow(s)); if("total" %in% names(s)) print(sum(s$total>=1000))'; sed -n 1,15p R/01_prepare_meta.R