# session cell 176 | 2026-11-02 11:11 | env: python | language: bash
# files written: []

F=/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/tool-results/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction-history/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction_0.txt
python3 - <<EOF
import re,json
s=open("$F").read()
i=s.find('AmphiBac_abundance_class_function.png')
# find the first python tool_use code that savefig this file
for m in re.finditer(r'savefig\(\\\\"AmphiBac_abundance_class_function.png', s):
    j=m.start(); st=s.rfind('"code":"',0,j); print(s[st:st+5000].replace('\\\\n','\n').replace('\\\\"','"')[:5000]); break
EOF