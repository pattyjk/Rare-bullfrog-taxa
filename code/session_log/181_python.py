# session cell 181 | 2026-11-02 11:12 | env: python | language: python
# files written: []

anyreads=X.sum(1)>0
S2=set(np.array(ids)[anyreads]); print(len(S2))
extra=S2-S
cnt=0
for k,s in fa.items():
    if I.loc[I.FastaID==k,'n_exact'].iat[0]>0: continue
    hit=False
    for q in (s,rc(s)):
        for i in range(len(q)-89):
            if q[i:i+90] in extra: hit=True;break
        if hit: break
    cnt+=hit
print(cnt)