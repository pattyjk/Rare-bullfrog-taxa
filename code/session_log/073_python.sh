# session cell 73 | 2026-10-17 16:56 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/16S_ribosomal_RNA.nog","sha256":"d43a675292dc0ab2ea3d7973b6e6eeb1330b5cfe74c9cf1d047bf2ccdff5e60f"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/16S_ribosomal_RNA.not","sha256":"1b9732e8c453f00b0443bc87308189c78b5caf3a7037b8754a4f9d25f98e5104"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/16S_ribosomal_RNA.nos","sha256":"a77067055ecf15c819a6963b75f35f63efd1e0e4dfe378cd8899f6effd1b54d2"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/taxonomy4blast.sqlite3","sha256":null,"oversize":true},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/16S_ribosomal_RNA.nto","sha256":"318081be517cace8aaee20d2914606848c743d51eed4a8e5c806940fdd6c3409"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/16S_ribosomal_RNA.ntf","sha256":"2c52678298b866ca5814764c7c288a46c7aa5ecfb69daaa24a7f1bbfd299ab20"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/16S_ribosomal_RNA.nhr","sha256":"59cc7930d36149ea8401d2b782d6f7f8c841b61bc21b487f44992281e2176ce7"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/taxdb.btd","sha256":null,"oversize":true},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/db/16S_ribosomal_RNA.tar.gz","sha256":null,"oversize":…[+4960 chars]

mkdir -p sal meta/ref meta/db
cp "<ARTIFACT asv_table.txt>" sal/asv_table.txt; cp "<ARTIFACT dna-sequences.fasta>" sal/asvs.fasta
head -2 sal/asv_table.txt | cut -c1-300; awk -F'\t' '{print NF}' sal/asv_table.txt | sort | uniq -c; grep -c ">" sal/asvs.fasta
grep -v ">" sal/asvs.fasta | awk '{print length($0)}' | sort | uniq -c | sort -k1nr | head -3
grep -v ">" sal/asvs.fasta | cut -c1-4 | sort | uniq -c | sort -k1nr | head -3
(cd meta/db && curl -s -m 600 -O https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz && curl -s -m 600 -O https://ftp.ncbi.nlm.nih.gov/blast/db/taxdb.tar.gz && for f in *.tar.gz; do tar xzf $f; done; ls | wc -l)
cd meta/ref; cat > fetch.py <<'E'
import urllib.request, urllib.parse, time, re
B="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for i in range(4):
        try: return urllib.request.urlopen(u,timeout=90).read().decode()
        except Exception: time.sleep(3)
    raise
for g in ["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas","Kineococcus","Rhodococcus","Paenarthrobacter","Pseudarthrobacter","Glutamicibacter"]:
    term=f'"{g}"[Organism] AND (16S ribosomal RNA[Title] OR 16S rRNA[Title]) AND 1000:1700[SLEN]'
    x=get(B+"esearch.fcgi?"+urllib.parse.urlencode({"db":"nuccore","term":term,"retmax":1500,"usehistory":"y"}))
    cnt=re.search(r"<Count>(\d+)",x).group(1); we=re.search(r"<WebEnv>(\S+?)</WebEnv>",x).group(1); qk=re.search(r"<QueryKey>(\d+)",x).group(1)
    n=min(int(cnt),1500); out=""
    for s in range(0,n,500):
        out+=get(B+"efetch.fcgi?"+urllib.parse.urlencode({"db":"nuccore","query_key":qk,"WebEnv":we,"rettype":"fasta","retmode":"text","retstart":s,"retmax":500})); time.sleep(0.4)
    open(f"{g}.fasta","w").write(out); print(g,cnt,out.count(">"),flush=True)
E
python fetch.py