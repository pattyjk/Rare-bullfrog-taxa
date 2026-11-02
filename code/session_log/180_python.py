# session cell 180 | 2026-11-02 11:12 | env: python | language: python
# files written: []

I['has_V4b']=I.has_V4|(I.n_exact>0)
det=(I.n_exact>0)|(I.n_near>0)
print(I.has_V4b.sum(), det.sum(), (~I.has_V4b).sum())
bd=meta.set_index('FastaID')
print([c for c in meta.columns if 'Bd' in c or 'Ref' in c or 'Growth' in c])