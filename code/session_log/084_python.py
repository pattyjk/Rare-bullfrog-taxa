# session cell 84 | 2026-10-17 17:50 | env: python | language: python
# files written: []


gen={"BF7":"24f22222-c59f-4660-a6c9-b1943330d6f0","BF3":"1d0f969c-960e-48fd-b8ab-e35e28955459","BF6":"cf52646c-30dc-48dd-baa0-2812086eee42","BF5":"75805d05-540f-4d38-b3b2-258490811841",
     "BF4":"ab02abc4-bcb5-4de5-8b24-d88541cb0e21","BF2":"61fdc6b1-6d5b-4a7b-844a-1d2d354aba93","BF1":"2faf35be-0a80-4e7b-9433-1cc8a4077080"}
rev=pr  # 515F
# also require region looks like 16S: check 806R downstream (GGACTACNVGGGTWTCTAAT) rc = ATTAGAWACCCBNGTAGTCC
r806=re.compile("ATTAGA[AT]ACCC[CGT][ACGT]GTAGTCC")
V4={}
for iso_,v in gen.items():
    seqs={};h=None
    for l in open(host.artifact_path(v)):
        l=l.strip()
        if l.startswith(">"): h=l[1:]; seqs[h]=[]
        elif h: seqs[h].append(l.upper())
    hits=[]
    for h,s in seqs.items():
        s="".join(s)
        for strand,qq in (("+",s),("-",s.translate(rc)[::-1])):
            for mm in pr.finditer(qq):
                seg=qq[mm.end():mm.end()+400]
                if r806.search(seg) and len(seg)>=120: hits.append((h[:30],len(s),seg[:120]))
    V4[iso_]=hits
    print(iso_, "16S V4 copies:",len(hits),"distinct 120bp:",len({x[2] for x in hits}), "contig lengths:",[x[1] for x in hits])
