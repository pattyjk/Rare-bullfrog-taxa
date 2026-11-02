# session cell 2 | 2026-10-15 19:40 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/all.txt","sha256":"633efafb0a4ba8e98685ab06dea53114fb485d572afd56ef39e0b75edfb04412"}]


full={k:"\n".join(t) for k,t in texts.items()}
with open("all.txt","w") as f:
    for k in ["HB58","BF246"]: f.write(f"#### {k}\n{full[k]}\n")
print(full["HB58"][3500:])
print("=================")
print(full["BF246"])
