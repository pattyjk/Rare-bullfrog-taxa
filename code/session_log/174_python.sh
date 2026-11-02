# session cell 174 | 2026-11-02 11:10 | env: python | language: bash
# files written: []

D=/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/tool-results/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction-history/893213d0-d199-4b99-9dfb-31fb80306f9a/
ls -la $D; grep -o "taxon[^\"]\{0,40\}\(majority\|any isolate\|frac_inhib\)[^\"]\{0,200\}" $D/*.txt | head -8; echo ---; grep -o "matched sOTU set[^\"]\{0,300\}" $D/*.txt | head -4; echo ---; grep -o "sOTU_set[^\"]\{0,200\}" $D/*.txt | head -5