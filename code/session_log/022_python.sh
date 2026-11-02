# session cell 22 | 2026-10-15 19:58 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/samples.txt","sha256":"e22a73a5006f8ea8b468a8546f91fbcbb75b15fafb5212b06c89636db3b793c2"}]

cd meta; python - <<'E'
import pandas as pd
m=pd.read_csv("map.txt",sep="\t",dtype=str)
print(m.shape, m["#SampleID"].nunique())
for c in ["host_order","captive_wild","external_habitat_type","lifestage","country","biogeo_realm","doi","target_subfragment","read_length_bp"]:
    print(c, m[c].value_counts().head(8).to_dict())
E
head -1 otu/amphib_only_otu_table.txt | tr '\t' '\n' > samples.txt; python - <<'E'
import pandas as pd
m=pd.read_csv("meta/map.txt" if False else "map.txt",sep="\t",dtype=str)
s=[l.strip() for l in open("samples.txt")]
print("otu samples",len(s),"in map",sum(x in set(m["#SampleID"]) for x in s))
E