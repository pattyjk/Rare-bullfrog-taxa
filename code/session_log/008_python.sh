# session cell 8 | 2026-10-15 19:46 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/map.txt","sha256":"9f4f7493f8fc47dbf37725d760fc1c1c26d66de86c897d44a680fdcef64056d5"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/otu.zip","sha256":"72f055a6a9e21f94ea207207505b2e234508c2b1d5233efab8f8c92050df3506"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/tax.txt","sha256":"bd9231ef1284c5eb702a5f41eb49a134dbf7893dd5c72270f63008add27d9332"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/otu/amphib_only_otu_table.txt","sha256":null,"oversize":true}]


T=$(python - <<'E'
E
)
cd "$(pwd)"; mkdir -p meta; cd meta
cp "<ARTIFACT meta_analysis_taxonomy.txt>" tax.txt
cp "<ARTIFACT amphib_only_map.txt>" map.txt
cp "<ARTIFACT amphib_only_otu_table.zip>" otu.zip
unzip -o -q otu.zip -d otu; find otu -type f | head; ls -la otu/* | head
wc -l tax.txt map.txt; head -3 tax.txt | cut -c1-300
head -2 map.txt | cut -c1-600
