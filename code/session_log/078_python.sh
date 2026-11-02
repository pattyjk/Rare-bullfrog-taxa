# session cell 78 | 2026-10-17 17:47 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amb/meta.tsv","sha256":"7ba987f352e3476ca3ce104f45a44d914f68d012a13000d2b9c32ed5cb58f994"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amb/seqs.fasta","sha256":"b274744fe67ee435034aed5100e04e993c91f7ccf7e9e9118f764e8c81a1729e"}]

mkdir -p amb; cp "<ARTIFACT AmphibBac_Metadata_2025.1-tested.tsv>" amb/meta.tsv; cp "<ARTIFACT AmphibBac_2025.1-tested.fasta>" amb/seqs.fasta
head -1 amb/meta.tsv | tr '\t' '\n' | nl | head -80; wc -l amb/meta.tsv; grep -c ">" amb/seqs.fasta; head -2 amb/seqs.fasta | cut -c1-200
python - <<'E'
import pandas as pd
m=pd.read_csv("amb/meta.tsv",sep="\t",dtype=str,low_memory=False)
print(m.shape)
for c in m.columns:
    if any(k in c.lower() for k in ["inhib","bd","result","function","facilit","genus","host","country"]):
        print(c, m[c].nunique(), m[c].value_counts().head(6).to_dict())
E