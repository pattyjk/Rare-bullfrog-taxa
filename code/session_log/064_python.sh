# session cell 64 | 2026-10-15 21:46 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/status.json","sha256":"47083fbc48d7547159862b19c6dba50f8b38698d54ed9e2746c0bc065b182b4c"}]

cd as; python - <<'E'
import json,time,urllib.request
jobs=[l.rstrip("\n").split("\t") for l in open("jobs.tsv")]
jobs={n:json.loads(j)["id"] for n,j in jobs}
t0=time.time(); st={}
while time.time()-t0<3000:
    for n,i in jobs.items():
        if st.get(n,{}).get("state") in ("done","failed"): continue
        try: st[n]=json.loads(urllib.request.urlopen(f"https://antismash.secondarymetabolites.org/api/v1.0/status/{i}",timeout=30).read())
        except Exception as e: pass
    if all(st.get(n,{}).get("state") in ("done","failed") for n in jobs): break
    time.sleep(60)
for n in jobs: s=st.get(n,{}); print(n,s.get("state"),s.get("status","")[:80],s.get("result_url",""))
json.dump(st,open("status.json","w"))
E