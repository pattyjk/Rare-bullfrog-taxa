# session cell 173 | 2026-11-02 11:10 | env: control-plane (host) | language: python
# files written: []

for q in ['taxon-level function majority', 'Ties -> most abundant sOTU', 'taxon inhibitory if', 'sOTU set', 'exact_sOTUs multiple']:
    r=host.archive.search(q,k=4)
    print('==',q)
    for h in (r if isinstance(r,list) else r.get('results',[])):
        print(str(h)[:600].replace('\\n',' '),'\n')