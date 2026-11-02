# session cell 16 | 2026-10-15 19:52 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/v4.py","sha256":"0fd73e3594e4a39ff7c3162e017ec283941ae8958e36ac4f6ea5c5250ffab730"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/v4refs.json","sha256":"59d6b3195e7cd5ceb93665bee016cdef093278745367f149dd3ac9e035073f51"}]


cd meta/ref
cat > v4.py <<'E'
import re, json
rc=str.maketrans("ACGTacgt","TGCAtgca")
pr=re.compile("GTG[CT]CAGC[AC]GCCGCGGTAA")
res={}
for g in ["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas"]:
    seqs={};h=None
    for l in open(f"{g}.fasta"):
        l=l.strip()
        if l.startswith(">"): h=l[1:]; seqs[h]=[]
        elif h: seqs[h].append(l.upper())
    v={}; nf=0
    for h,s in seqs.items():
        s="".join(s)
        for q in (s, s.translate(rc)[::-1]):
            m=pr.search(q)
            if m and len(q)>=m.end()+90:
                frag=q[m.end():m.end()+90]
                if set(frag)<=set("ACGT"): v[h]=frag
                break
        else: nf+=1
    res[g]=v; print(g,len(seqs),"V4 extracted",len(v),"unique",len(set(v.values())))
json.dump(res,open("v4refs.json","w"))
E
python v4.py
