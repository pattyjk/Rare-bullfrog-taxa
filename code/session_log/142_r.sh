# session cell 142 | 2026-10-18 11:48 | env: r | language: bash
# files written: []

cd amphib_R; Rscript -e 'library(data.table); s<-fread("output/meta_genus_summary.csv"); options(width=200); print(dcast(s,genus~method,value.var="prevalence_pct"),digits=3); print(dcast(s,genus~method,value.var="mean_rel_pct"),digits=3)'