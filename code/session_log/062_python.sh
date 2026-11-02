# session cell 62 | 2026-10-15 21:33 | env: python | language: bash
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/as/jobs.tsv","sha256":"1abb119a13ec0dca81bcf188589839f63205f04e8730b3af337f900f414d7e2a"}]

cd as; : > jobs.tsv
for f in *.min1k.fasta; do
  n=${f%.min1k.fasta}
  r=$(curl -s -m 300 -X POST https://antismash.secondarymetabolites.org/api/v1.0/submit \
     -F "seq=@$f" -F knownclusterblast=true -F cc_mibig=true -F clusterhmmer=true -F asf=true -F rre=true -F tta=true -F pfam2go=true \
     -F hmmdetection_strictness=relaxed -F genefinder=prodigal)
  echo -e "$n\t$r" | tee -a jobs.tsv
  sleep 2
done