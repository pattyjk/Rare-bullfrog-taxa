# session cell 108 | 2026-10-17 21:32 | env: python | language: python
# files written: []


mpg["hg"]=mpg.host_genus.str.lower().str.strip(); mpg["hf"]=mpg.host_family.str.lower().str.strip()
g2f=mpg.groupby("hg").hf.agg(lambda x:x.value_counts().index[0]).to_dict()
print("meta families:",sorted(mpg.hf.unique()))
m["hg"]=m.hsp.str.split().str[0]
mg=m[m.sOTUs!=""]
miss_g=mg[~mg.hg.isin(g2f)].hg.value_counts(); print("genera not in meta:",miss_g.to_dict())
