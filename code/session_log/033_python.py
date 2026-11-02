# session cell 33 | 2026-10-15 20:52 | env: python | language: python
# files written: []


import collections
ids=[s for s in exm[exm.ref_genera.str.contains("Rhodococcoides")].sOTU]
for s in ids:
    out={}
    for g in ["Rhodococcoides","Rhodococcus"]:
        c=collections.Counter(" ".join(h.split()[1:3]) for h,v in refs[g].items() if v==s)
        out[g]=dict(c.most_common(6))
    print(s in present, int(cnt.loc[s,good].sum()) if s in present else 0, out)
