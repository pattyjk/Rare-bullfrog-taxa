# session cell 175 | 2026-11-02 11:11 | env: python | language: bash
# files written: []

F=/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/tool-results/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction-history/893213d0-d199-4b99-9dfb-31fb80306f9a/compaction_0.txt
python3 - <<EOF
import re
s=open("$F").read()
for pat in ['sOTU_set','taxon_func','majority','any_inhib','frac_inhib >= 0.5','pct_inhib']:
    idx=[m.start() for m in re.finditer(re.escape(pat),s)]
    print('##',pat,len(idx))
    for i in idx[:3]: print(s[max(0,i-300):i+300].replace('\\\\n','\n')[:600]); print('....')
EOF