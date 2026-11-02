# session cell 37 | 2026-10-15 20:58 | env: python | language: python
# files written: []


rc=str.maketrans("ACGT","TGCA")
gens=["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas","Rhodococcus","Paenarthrobacter","Pseudarthrobacter","Glutamicibacter"]
ref120={}  # variant -> {genus: Counter(species)}
nref={}
for g in gens:
    seqs={};h=None
    for l in open(f"meta/ref/{g}.fasta"):
        l=l.strip()
        if l.startswith(">"): h=l[1:]; seqs[h]=[]
        elif h: seqs[h].append(l.upper())
    n=0
    for h,s in seqs.items():
        s="".join(s)
        for qq in (s,s.translate(rc)[::-1]):
            mm=pr.search(qq)
            if mm and len(qq)>=mm.end()+120:
                fr=qq[mm.end():mm.end()+120]
                if set(fr)<=set("ACGT"):
                    sp=" ".join(h.split()[1:3]); ref120.setdefault(fr,{}).setdefault(g,collections.Counter())[sp]+=1; n+=1
                break
    nref[g]=n
print(nref, "unique 120bp variants:",len(ref120))
V=list(ref120); Va=np.frombuffer("".join(V).encode(),dtype=np.uint8).reshape(len(V),120)
rows=[]
for aid,s in q.items():
    a=np.frombuffer(s.encode(),dtype=np.uint8)
    d=(Va!=a).sum(1); dmin=d.min()
    if dmin<=2:
        hits=[V[i] for i in np.where(d==dmin)[0]]
        gset={}; 
        for v in hits:
            for g,c in ref120[v].items(): gset.setdefault(g,collections.Counter()).update(c)
        rows.append(dict(ASV_ID=aid,mismatches=int(dmin),pct_identity=round(100*(120-dmin)/120,2),
            matched_genera=";".join(sorted(gset)),
            top_ref_species="; ".join(f"{sp} ({n})" for g in sorted(gset) for sp,n in gset[g].most_common(3)),
            sequence=s))
A=pd.DataFrame(rows)
print(A.groupby("mismatches").size().to_dict())
print(A[A.mismatches==0].matched_genera.value_counts().to_dict())
