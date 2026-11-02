# session cell 167 | 2026-10-28 18:35 | env: python | language: python
# files written: []

r=host.artifacts(limit=500)
arts=r['artifacts'] if isinstance(r,dict) else r
import collections
rows=sorted([(a['filename'],a.get('latest_version_id'),a.get('created_at','')[:10]) for a in arts])
print(len(rows))
for x in rows: print(*x, sep=' | ')