# session cell 58 | 2026-10-15 21:29 | env: python | language: bash
# files written: []

which docker podman brew 2>&1; docker info 2>&1 | head -3; curl -s -m 15 -o /dev/null -w "%{http_code}\n" https://antismash.secondarymetabolites.org/api/v1.0/version; curl -s -m 15 -o /dev/null -w "%{http_code}\n" https://pypi.org/pypi/antismash/json