# session cell 178 | 2026-11-02 11:12 | env: python | language: python
# files written: []

t0=time.time()
def read_fasta(p):
    d={};k=None;b=[]
    for l in open(p):
        l=l.strip()
        if l.startswith('>'):
            if k: d[k]=''.join(b)
            k=l[1:].split()[0];b=[]
        else: b.append(l.upper())
    if k: d[k]=''.join(b)
    return d
fa=read_fasta(host.artifact_path('a7223533-0656-479a-a562-6f2627103fea'))
meta=pd.read_csv(host.artifact_path('e9e6f252-b2f8-43f7-8080-ccdeecc5c8ca'),sep='\t',low_memory=False)
comp=str.maketrans('ACGTRYKMBDHVN','TGCAYRMKVHDBN')
rc=lambda s:s.translate(comp)[::-1]
P=re.compile('GTG[CT]CAGC[AC]GCCGCGGTAA')
def win(s,n=90):
    for q in (s,rc(s)):
        m=P.search(q)
        if m and m.end()+n<=len(q): return q[m.end():m.end()+n]
    return None
S=set(st.sOTU); ab=dict(zip(st.sOTU,st.mean_rel_pct))
M=np.frombuffer(''.join(st.sOTU).encode(),dtype=np.uint8).reshape(-1,90)
rec=[]
for k,s in fa.items():
    ex=set()
    for q in (s,rc(s)):
        for i in range(len(q)-89):
            w=q[i:i+90]
            if w in S: ex.add(w)
    w=win(s); dmin=None; near=[]
    if not ex and w is not None and re.fullmatch('[ACGTN]+',w):
        d=(M!=np.frombuffer(w.encode(),dtype=np.uint8)).sum(1); dmin=int(d.min())
        if dmin<=3: near=list(st.sOTU.values[d==dmin])
    rec.append((k,w is not None,';'.join(sorted(ex)),len(ex),dmin,';'.join(near),len(near)))
I=pd.DataFrame(rec,columns=['FastaID','has_V4','exact_sOTUs','n_exact','nearest_mm','near_sOTUs','n_near'])
print(len(I), I.has_V4.sum(), ((I.n_exact>0)|(I.n_near>0)).sum(), ((I.n_exact>0)|(I.n_near>0))[I.has_V4].sum(), round(time.time()-t0))