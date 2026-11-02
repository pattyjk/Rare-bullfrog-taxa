# session cell 60 | 2026-10-15 21:32 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF6.fasta","sha256":"a90a239dddd2096fca9263b3e9534e279c903555ccc4bf031ca7c2432198c641"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF7.fasta","sha256":"1bf833baa0216845cd8e479206b00d3a642fb8bf3b53bbcf822ededdd4b798cc"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF4.fasta","sha256":"a4b79b91d74a9542a50b9805e468199b5a35f4c42e6ccd2cca58dbb8b23e8a31"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF1.fasta","sha256":"f39a6c325dede313c69fbb9a1f3dcc23a7482dfa465cdadaf06e5a42152cd706"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF5.fasta","sha256":"413889085387ea022811e08a328f71b59ab76f49c934e282ab5cbd850f496352"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF3.fasta","sha256":"933f4b11fcf651976fc0441c5478f5e9cf405a1a16e2e6cdab350a62861a9417"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/BF2.fasta","sha256":"a412afa41b8cd036264e28f8121fc26524a4102b3ed08161deecfd7f3ec9dfad"}]

mkdir -p as && cd as
cp "<ARTIFACT BF1.fasta>" BF1.fasta
cp "<ARTIFACT BF2.fasta>" BF2.fasta
cp "<ARTIFACT BF3.fasta>" BF3.fasta
cp "<ARTIFACT BF4.fasta>" BF4.fasta
cp "<ARTIFACT BF5.fasta>" BF5.fasta
cp "<ARTIFACT BF6.fasta>" BF6.fasta
cp "<ARTIFACT BF7.fasta>" BF7.fasta
for f in *.fasta; do echo "$f $(grep -c '>' $f) contigs, $(grep -v '>' $f | tr -d '\n' | wc -c) bp, hdr: $(head -1 $f | cut -c1-50)"; done