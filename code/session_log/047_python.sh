# session cell 47 | 2026-10-15 21:10 | env: python | language: bash
# files written: []

f="<ARTIFACT asv_table.txt>"; head -3 "$f" | cut -c1-400; wc -l "$f"; awk -F'\t' '{print NF}' "$f" | sort | uniq -c | head