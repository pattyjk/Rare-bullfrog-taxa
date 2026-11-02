# session cell 157 | 2026-10-18 12:51 | env: r | language: bash
# files written: []

cd amphib_R; ls -la /tmp/fetchtest; grep -c ">" /tmp/fetchtest/*.fasta 2>/dev/null; grep -c ">" data/raw/ncbi_refs/Psychromicrobium.fasta
cat output/Bd_variance_partitioning.csv
Rscript -e 'library(data.table); d<-fread("output/sOTU_Bd_function_variance.csv"); print(names(d)); print(nrow(d)); s<-grep("sd",names(d),value=TRUE,ignore.case=TRUE); print(s); print(summary(d[[s[1]]]))'