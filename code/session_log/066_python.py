# session cell 66 | 2026-10-15 21:58 | env: python | language: python
# files written: []


import json,glob,os
J=json.load(open("as/res/BF1/BF1.min1k.json"))
print(J.keys(), J.get("version"))
rec=[r for r in J["records"] if r.get("areas")][0]
print(rec.keys()); print(rec["areas"][0])
mods=rec["modules"]; print(list(mods.keys()))
kc=mods.get("antismash.modules.clusterblast",{}); print(kc.keys())
k=kc["knowncluster"]; print(k.keys()); r0=k["results"][0]; print(r0.keys()); print(json.dumps(r0["ranking"][0])[:1500] if r0["ranking"] else r0)
