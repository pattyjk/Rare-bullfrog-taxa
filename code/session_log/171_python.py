# session cell 171 | 2026-11-02 11:10 | env: python | language: python
# files written: []

r=host.artifacts(limit=500); arts=r['artifacts'] if isinstance(r,dict) else r
for a in sorted(arts,key=lambda a:a.get('created_at',''))[-40:]:
    if a.get('root_frame_id')!='893213d0-d199-4b99-9dfb-31fb80306f9a': print(a['created_at'][:16],a['filename'],a['latest_version_id'],a.get('root_frame_id','')[:8])