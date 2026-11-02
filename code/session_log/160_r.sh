# session cell 160 | 2026-10-18 12:52 | env: r | language: bash
# files written: []

cd amphib_R; Rscript -e '
library(data.table); S<-fread("output/sOTU_Bd_function_variance.csv"); V<-S[n_scored>=5]
cat(nrow(V), median(V$sd_score), 100*mean(V$max_score>50 & V$min_score<0), "\n")
print(V[, .(n=.N, sd=median(sd_score)), by=abundance_class]); print(kruskal.test(sd_score~factor(abundance_class), V)$p.value)
A<-S[n_isolates>=3]; cat(nrow(A), median(A$label_agreement_pct), 100*mean(A$label_agreement_pct==100), "\n")'
grep -o '"Bd_function[^)]*' -m3 /Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/tool-results/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction-history/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction_0.txt | head -3; grep -n "175 sOTUs\|n = 175\|175)" /Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/tool-results/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction-history/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction_0.txt | head -5