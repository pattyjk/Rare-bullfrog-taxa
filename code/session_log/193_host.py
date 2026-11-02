# session cell 193 | 2026-11-02 11:38 | env: control-plane (host) | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/handoff/exec_log.json","sha256":"f0b7bf280fc243040d120accc83eb6023fa2e204839d80b7acee92b0e0d35f3d"}]

import json, os
F='893213d0-d199-4b99-9dfb-31fb80306f9a'
r=host.query(f"SELECT count(*) FROM execution_log WHERE frame_id='{F}'"); print(r['rows'])
rows=[]; off=0
while True:
    q=host.query(f"SELECT cell_index, created_at, conda_env, language, exit_status, files_written, source FROM execution_log WHERE frame_id='{F}' ORDER BY created_at LIMIT 200 OFFSET {off}")
    rows+=q['rows']; 
    if len(q['rows'])<200: break
    off+=200
os.makedirs('handoff',exist_ok=True)
json.dump({'columns':['cell_index','created_at','conda_env','language','exit_status','files_written','source'],'rows':rows},open('handoff/exec_log.json','w'))
from collections import Counter
print(len(rows), Counter((r[3],r[4]) for r in rows))