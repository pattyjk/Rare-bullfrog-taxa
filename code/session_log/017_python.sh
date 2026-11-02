# session cell 17 | 2026-10-15 19:53 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/v4b.py","sha256":"a00d9e3567743032ec72a024f4f8a6808f47f33b75e5081b0bd881ae20a98330"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Rhodococcus.fasta","sha256":"2844057ce2f2fb700cdb392f21064f60c1b24fc9a020bbb3c95a71880143b420"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/fetch2.py","sha256":"62608e021779d958054fda49365ad07a4fa0868c82f7628863725f13a0fdec65"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/v4refs.json","sha256":"fbe514c5a495d7de8ae508d6a55a9148fecc899217bdec7cbb23806aaf55847b"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Paenarthrobacter.fasta","sha256":"460324d1918d42cda587980b17ef8939a591b753fa977792b263ea4f87ba69ee"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Pseudarthrobacter.fasta","sha256":"a46c79e4d3c25fd48f72c78ca28b659008a132c557548bf9e9a7c95be2891870"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/ref/Glutamicibacter.fasta","sha256":"2e092f9f2362dd55c1d12ee48eb28e2c4cb728bb66fcf1f87a31563534205ab7"}]


cd meta/ref
sed -e 's/for g in \["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas"\]/for g in ["Rhodococcus","Paenarthrobacter","Pseudarthrobacter","Glutamicibacter"]/' fetch.py > fetch2.py
python fetch2.py
sed -e 's/for g in \["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas"\]/for g in ["Rhodococcoides","Psychromicrobium","Arthrobacter","Cryobacterium","Ancylobacter","Sphingomonas","Rhodococcus","Paenarthrobacter","Pseudarthrobacter","Glutamicibacter"]/' v4.py > v4b.py
python v4b.py
python - <<'E'
import json,collections
r=json.load(open("v4refs.json"))
U={g:set(v.values()) for g,v in r.items()}
gs=list(U)
for i,a in enumerate(gs):
    for b in gs[i+1:]:
        o=len(U[a]&U[b])
        if o: print(a,b,"shared V4 variants:",o)
E
