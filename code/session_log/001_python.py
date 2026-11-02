# session cell 1 | 2026-10-15 19:40 | env: python | language: python
# files written: []


import fitz
texts={}
for k,v in vids.items():
    d=fitz.open(host.artifact_path(v)); texts[k]=[pg.get_text() for pg in d]
for k,t in texts.items(): print(k,len(t),sum(len(x) for x in t))
print(texts["HB58"][0][:4000])
