# session cell 146 | 2026-10-18 12:07 | env: r | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/figures/sOTU_Bd_variance.png","sha256":"d7656999d076a1e92a9c67be9400776edec88c6786a602e7cd4ecd7520782865","preview_sha256":"49f8ea743934fa45972493b4a30f4feba103415d557f480cc28b0bef933238d3","preview_content_type":"image/jpeg"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/sOTU_Bd_function_variance.csv","sha256":"4e0cf702cf9b535808f31112957bc3daff679135a95687c6e584cdda4167e5be"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/output/Bd_variance_partitioning.csv","sha256":"2394ad8f158ca9c31e18e97d2e683953003bdf0da287f79b33034d2a37443681"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/amphib_R/R/05_sOTU_function_variance.R","sha256":"08ac4966263b3c421a0b77c859eb849cba1aaacae7f57168a69e02e25afc1e1c"}]

cd amphib_R
python3 - <<'E'
p="R/05_sOTU_function_variance.R"; s=open(p).read()
s=s.replace('title = "Spread is similar across classes"','title = "Similar spread by class"').replace('title = "Label agreement within sOTUs"','title = "Label agreement"')
s=s.replace('labels = c("Dominant\\n>=1%", "Interm.\\n0.1-1%", "Rare\\n<0.1%")','labels = c("Dominant\\n\\u22651%", "Interm.\\n0.1\\u20131%", "Rare\\n<0.1%")')
s=s.replace('colc <- c("Dominant (>=1%)" = "#08306B", "Intermediate (0.1-1%)" = "#4292C6", "Rare (<0.1%)" = "#BDBDBD")','colc <- c("Dominant (>=1%)" = "#08306B", "Intermediate (0.1-1%)" = "#4292C6", "Rare (<0.1%)" = "#BDBDBD")\nclab <- c("Dominant (\\u22651%)", "Intermediate (0.1\\u20131%)", "Rare (<0.1%)")')
s=s.replace('scale_colour_manual(values = colc)','scale_colour_manual(values = colc, labels = clab)')
s=s.replace('"sOTUs with >=5 scored isolates (n = %d), ranked by median"','"sOTUs with \\u22655 scored isolates (n = %d), ranked by median"')
open(p,"w").write(s)
E
Rscript R/05_sOTU_function_variance.R >/dev/null 2>&1; ls -la figures output | awk '{print $5,$9}'