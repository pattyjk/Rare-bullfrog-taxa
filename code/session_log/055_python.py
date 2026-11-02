# session cell 55 | 2026-10-15 21:27 | env: python | language: python
# files written: []


import hashlib
new={"BF1":"97b1b65c-daa4-4850-86b8-2f38fc719071","BF2":"44fa9764-1096-4983-add6-8291ffac2a31","BF3":"584106f1-c286-4c9d-9bad-8e5c29d470ee","BF5":"17ed7aff-a1e2-4a24-bd3c-d5eeff718c36","BF4":"2b208ffc-3d6f-443c-a1a5-a50eecc8135e","BF6":"828b5838-8ef7-4b48-a223-8139d5eea73e","BF7":"8cc80d95-cf2f-4e00-b092-73e0d64a87ce"}
def md5(v): return hashlib.md5(open(host.artifact_path(v),"rb").read()).hexdigest()
for k,v in new.items():
    print(k, "identical to earlier upload" if k in vids and md5(v)==md5(vids[k]) else "NEW/changed")
fullN={k:"\n".join(pg.get_text() for pg in fitz.open(host.artifact_path(v))) for k,v in new.items()}
print(fullN["BF3"][:3500])
