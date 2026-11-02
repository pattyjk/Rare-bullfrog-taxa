# session cell 7 | 2026-10-15 19:41 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/MiGA_isolate_summary.csv","sha256":"c82b3b777831349dd43a0ce7468ed7b8a711bf65262ac8db9f96ffe6b5a9259a"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/MiGA_isolate_summary.xlsx","sha256":"7dcf591c0f13eccc01e9006ec5e0217ff2808e2c081b74baaa7205fd8ceea35c"}]


print(df.loc[df.file=="BF4old","novelty"].iloc[0])
df.loc[df.file=="HB58","interpretation"]="Acinetobacter oleivorans; borderline (ANI 95.25%, just above the 95% species cutoff). Folder label 'Kineo' does not match."
out=df[cols]; out.to_csv("MiGA_isolate_summary.csv",index=False); out.to_excel("MiGA_isolate_summary.xlsx",index=False)
