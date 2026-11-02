# session cell 18 | 2026-10-15 19:54 | env: python | language: bash
# files written: []

cd meta/otu; awk -F'\t' '{print NF}' amphib_only_otu_table.txt | uniq -c | head -20; awk -F'\t' 'NF!=2437{print NR": "substr($0,1,200)}' amphib_only_otu_table.txt | head -5; wc -l amphib_only_otu_table.txt