# session cell 139 | 2026-10-18 11:46 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/data/derived/sample_totals.csv","sha256":"914ee86a973e8d2fc2072e7e830cc6499a02cdf710ddac63f1c6c9f979bf7d06"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/data/derived/sOTU_stats.csv","sha256":"e320371b80626ab8644396681fe3e0dae4d892fefab86054d67470e64e617c6e"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/data/derived/rel_abund_sparse.rds","sha256":"607ce4ed8dbffdbac7781d8a7d08bbf174ce0f79d00cbf1c9e309f6444781472"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/01_prepare_meta.R","sha256":"e1ac57e1d086ce01128d2cd623f3671c2e78145a64b2cdca66c0a54b45f0526c"}]

cd amphib_R
python3 - <<'E'
p="R/01_prepare_meta.R"; s=open(p).read()
s=s.replace('st <- data.table(sOTU = rownames(R), reads = rowSums(M[rownames(R), keep]),\n                 n_samples = diff(R@p)[0] , stringsAsFactors = FALSE)','st <- data.table(sOTU = rownames(R), reads = rowSums(M[rownames(R), keep]))')
s=s.replace('st[, max_rel_pct := apply(R, 1, max)]','mx <- numeric(nrow(R)); tm <- tapply(R@x, R@i + 1, max); mx[as.integer(names(tm))] <- tm\nst[, max_rel_pct := mx]')
open(p,"w").write(s)
E
time Rscript R/01_prepare_meta.R 2>&1 | tail -5
ls -la data/derived | awk '{print $5,$9}'
head -3 data/derived/sOTU_stats.csv