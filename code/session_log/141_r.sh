# session cell 141 | 2026-10-18 11:47 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/data/derived/ncbi_ref_v4.csv","sha256":"bbbff8e210a8a651f81a2b2f2d39c09ed22731cd251b5c5e5e5cce3483a0ba96"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/meta_genus_by_host_metadata.csv","sha256":"5e147d5ee473c048850a2950f096aeada330acef7f3cb8fb49b495fb34c83e2f"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/meta_genus_relabund_by_sample.csv","sha256":"cf7913d853299d9bb4f8b4e8590798c384c83847224d4b4d8a1fe35a1a7371cf"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/meta_genus_summary.csv","sha256":"aeaf522a574e9688bd4413159c126d5de589ae86265100ba4f5ccd7be4a23c4c"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/meta_target_sOTUs.csv","sha256":"d467c8a33f6ec3e12253e3f4f6d6aea56156079a919629869256dc73f53e257e"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/02_meta_genus_abundance.R","sha256":"267a65d24b39ce5087019e59437c8b6c3218d5cab02c2c2aaaccc6fdbb0a1a97"}]

cd amphib_R
python3 - <<'E'
p="R/02_meta_genus_abundance.R"; s=open(p).read()
s=s.replace("ex <- refs[!is.na(v4_90), .(n_refs = .N, genera = paste(sort(unique(lab(ref_genus))), collapse = \";\")), by = .(v4_90, g = lab(ref_genus))]",
 "ex <- refs[!is.na(v4_90), .(n_refs = .N), by = .(v4_90, tg = lab(ref_genus))]")
s=s.replace("    e <- ex[g == ..g]","    e <- ex[tg == g]")
open(p,"w").write(s)
E
time Rscript R/02_meta_genus_abundance.R 2>&1 | tail -15