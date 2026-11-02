# session cell 179 | 2026-11-02 11:12 | env: python | language: python
# files written: []

from collections import Counter
c=Counter(''.join(fa.values())); print({k:v for k,v in c.items() if k not in 'ACGT'})