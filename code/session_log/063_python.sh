# session cell 63 | 2026-10-15 21:34 | env: python | language: bash
# files written: []

cd as; sleep 60; while IFS=$'\t' read n j; do id=$(echo $j | python -c "import sys,json;print(json.load(sys.stdin)['id'])"); echo "$n $(curl -s -m 30 https://antismash.secondarymetabolites.org/api/v1.0/status/$id)"; done < jobs.tsv