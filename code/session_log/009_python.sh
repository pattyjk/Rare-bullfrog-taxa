# session cell 9 | 2026-10-15 19:47 | env: python | language: bash
# files written: []


cd meta
head -c 600 otu/amphib_only_otu_table.txt; echo; echo ---
sed -n 2p otu/amphib_only_otu_table.txt | cut -c1-400; echo ---
awk -F'\t' 'NR<=3{print NF}' otu/amphib_only_otu_table.txt
for g in Rhodococcoides Rhodococcus Psychromicrobium Arthrobacter Cryobacterium Ancylobacter Sphingomonas Microterricola Starkeya; do echo "$g $(cut -f6 tax.txt | grep -c "; $g;\|; $g$")"; done
cut -f6 tax.txt | grep -m3 "Rhodococc"
cut -f6 tax.txt | awk -F'; ' '{print NF}' | sort | uniq -c | sort -k1nr | head
