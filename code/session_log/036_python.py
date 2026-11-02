# session cell 36 | 2026-10-15 20:58 | env: python | language: python
# files written: []


p=host.artifact_path("5d7e73fe-7eea-4eb8-af90-cffc1d6c04b7")
q={};h=None
for l in open(p):
    l=l.strip()
    if l.startswith(">"): h=l[1:]; q[h]=[]
    elif h: q[h].append(l.upper())
q={k:"".join(v) for k,v in q.items()}
L=pd.Series({k:len(v) for k,v in q.items()})
print(len(q), L.describe().round(0).to_dict(), L.value_counts().head(5).to_dict())
print(list(q.items())[:2])
pr=re.compile("GTG[CT]CAGC[AC]GCCGCGGTAA")
print("with 515F primer:",sum(bool(pr.search(s)) for s in q.values()), "start TAC:",sum(s.startswith("TAC") for s in q.values()))
