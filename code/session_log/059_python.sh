# session cell 59 | 2026-10-15 21:31 | env: python | language: bash
# files written: []

curl -s -m 20 https://antismash.secondarymetabolites.org/api/v1.0/version; echo; curl -s -m 20 https://antismash.secondarymetabolites.org/api/v1.0/stats | head -c 400; echo; curl -s -m 20 https://antismash.secondarymetabolites.org/ | grep -oiE "antismash [0-9]+\.[0-9]+[^<]{0,20}" | head -3