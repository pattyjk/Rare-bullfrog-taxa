# session cell 15 | 2026-10-15 19:52 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/fetch.py","sha256":"0fc1d3a2d61fda29224ad467fd5868a43def766210a202f8171b8bd43875f24e"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Sphingomonas.fasta","sha256":"993dbc7e44154ff811b26840efd98e5fe6ab0ab2a01dc836985993ca8f03b0fd"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Psychromicrobium.fasta","sha256":"e18691ef10746e56011b62f2916da0a932d4b400e2488cf421d1a3d9d92b98a9"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Ancylobacter.fasta","sha256":"bb9508b32216a154c3e0714b3a515b3746d68fb51bbc6aa51780b22e9856bd0e"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Cryobacterium.fasta","sha256":"71ad08a8e225d826338329b9304994ea53d7fa6c0364a9834b0bd8b71b3e5eda"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Arthrobacter.fasta","sha256":"a3a6bef6404bd3049e00367101c9e399ee1dbde80486c38bb1bd7c6f6e6d4362"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Rhodococcoides.fasta","sha256":"b037edc6a554b73dbba97f48f047cdd8d16a47460d56bcfb3316b3c8a097519d"}]


mkdir -p meta/ref && cd meta/ref
cat > fetch.py <<'E'
import urllib.request, urllib.parse, time, re, sys
B="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for i in range(3):
        try: return urllib.request.urlopen(u,timeout=60).read().decode()
        except Exception as e: time.sleep(2)
    raise
for g in ["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas"]:
    term=f'"{g}"[Organism] AND (16S ribosomal RNA[Title] OR 16S rRNA[Title]) AND 1000:1700[SLEN]'
    x=get(B+"esearch.fcgi?"+urllib.parse.urlencode({"db":"nuccore","term":term,"retmax":1500,"usehistory":"y"}))
    cnt=re.search(r"<Count>(\d+)",x).group(1); we=re.search(r"<WebEnv>(\S+?)</WebEnv>",x).group(1); qk=re.search(r"<QueryKey>(\d+)",x).group(1)
    n=min(int(cnt),1500); out=""
    for s in range(0,n,500):
        out+=get(B+"efetch.fcgi?"+urllib.parse.urlencode({"db":"nuccore","query_key":qk,"WebEnv":we,"rettype":"fasta","retmode":"text","retstart":s,"retmax":500})); time.sleep(0.4)
    open(f"{g}.fasta","w").write(out); print(g,cnt,out.count(">"),flush=True)
E
python fetch.py
