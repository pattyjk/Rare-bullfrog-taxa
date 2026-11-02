# session cell 86 | 2026-10-17 17:51 | env: python | language: python
# files written: []


seqs={};h=None
for l in open(host.artifact_path(gen["BF1"])):
    l=l.strip()
    if l.startswith(">"): h=l[1:]; seqs[h]=[]
    else: seqs[h].append(l)
cov=pd.Series({h:float(re.search(r"cov_([\d.]+)",h).group(1)) for h in seqs}); ln=pd.Series({h:sum(map(len,v)) for h,v in seqs.items()})
big=ln>=1000
print("BF1 contigs>=1kb weighted median cov:", round(np.average(cov[big],weights=ln[big]),1), "| hit contigs:",[(x[0],) for x in V4["BF1"]])
for x in V4["BF1"]:
    hh=[h for h in seqs if h.startswith(x[0])][0]; print(hh, ln[hh])
# coverage distribution: fraction of >=1kb assembly with cov < 10
print("bp in contigs cov<10:", int(ln[big&(cov<10)].sum()), "of", int(ln[big].sum()))
print(cov[big].describe().round(1).to_dict())
