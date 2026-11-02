# session cell 182 | 2026-11-02 11:13 | env: python | language: python
# files written: []

amb=I[(I.n_exact==0)&I.has_V4].copy()
amb['w']=[win(fa[k]) for k in amb.FastaID]
amb=amb[~amb.w.str.fullmatch('[ACGTN]+')]
d3=[int((M!=np.frombuffer(w.encode(),dtype=np.uint8)).sum(1).min()) for w in amb.w]
print(len(amb), sum(x<=3 for x in d3))