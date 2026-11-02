# session cell 25 | 2026-10-15 19:59 | env: python | language: python
# files written: [{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/totals.pkl","sha256":"c036e53d75c8fda4955e52d210e7484f60a8733e88ba3298b121d130900ef63e"},{"path":"/Users/p.kearns/.claude-science/orgs/02812373-138d-4a58-874b-33e46fa50871/workspaces/893213d0-d199-4b99-9dfb-31fb80306f9a/meta/sub_counts.pkl","sha256":"93b5bce6e553f99a4fcdf1423451397bbc5e1599b709e842d7adfd00557ad077"}]


sub.to_pickle("meta/sub_counts.pkl"); totals.to_pickle("meta/totals.pkl")
print("seq rows",nseq,"numeric rows",nnum,"bad",len(bad),"sub",sub.shape)
print("bad in keep:",sum(b[0] in keep for b in bad), "bad ids len90:",sum(len(b[0])==90 for b in bad), [ (b[0][:20],b[1]) for b in bad[:5]])
print({k:float(v) for k,v in totals.sum().items()})
T=totals.seq_reads+totals.numericID_reads
print((totals.numericID_reads/T).describe().round(3).to_dict())
print("samples with numeric-ID reads>50%:",(totals.numericID_reads/T>0.5).sum())
print(T.describe().round(0).to_dict(), (T<1000).sum())
