# session cell 20 | 2026-10-15 19:58 | env: python | language: bash
# files written: []

cd meta/otu; f=amphib_only_otu_table.txt
awk -F'\t' 'NR>1 && length($1)!=90{print NR, NF, substr($0,1,120)}' $f | head -8
awk -F'\t' 'NR>1 && length($1)!=90{print NR}' $f | awk 'NR==1{a=$1} END{print "first",a,"last",$1}'
awk -F'\t' 'NR>1 && length($1)==90{print NR}' $f | tail -1
awk -F'\t' 'NR>1 && length($1)!=90 && NF!=2437' $f | wc -l