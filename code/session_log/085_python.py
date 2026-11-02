# session cell 85 | 2026-10-17 17:50 | env: python | language: python
# files written: []


# bullfrog & salamander
def load_asv(tv,fv,skip=0):
    t=pd.read_csv(host.artifact_path(tv),sep="\t",index_col=0,skiprows=skip); qq={};h=None
    for l in open(host.artifact_path(fv)):
        l=l.strip()
        if l.startswith(">"): h=l[1:]; qq[h]=""
        else: qq[h]+=l.upper()
    return t,qq
bft,bfq=load_asv("b3169217-410e-4dad-bff1-4c925a06efb7","5d7e73fe-7eea-4eb8-af90-cffc1d6c04b7")
sat,saq=tab,q
def asv_stats(t,qq,seq120):
    T_=t.sum(); incl=T_>=1000; R_=t.loc[:,incl].div(T_[incl],axis=1)*100
    s2id={v:k for k,v in qq.items()}
    ids=[k for k,v in qq.items() if v==seq120]
    if ids:
        r_=R_.loc[ids].sum(); return dict(match="exact 120 bp",ASV=ids[0][:10],prev=round(100*(r_>0).mean(),1),mean=round(r_.mean(),4),max=round(r_.max(),3),mm=0)
    # nearest
    A=np.frombuffer("".join(qq.values()).encode(),dtype=np.uint8).reshape(len(qq),120); a_=np.frombuffer(seq120.encode(),dtype=np.uint8)
    d=(A!=a_).sum(1); i=d.argmin(); k=list(qq)[i]; r_=R_.loc[k]
    return dict(match=f"nearest ({d[i]} mismatch{'es' if d[i]>1 else ''})",ASV=k[:10],prev=round(100*(r_>0).mean(),1),mean=round(r_.mean(),4),max=round(r_.max(),3),mm=int(d[i]))
out=[]
for iso_,hits in V4.items():
    for ci,(ctg,clen,s120) in enumerate(hits,1):
        s90=s120[:90]
        r={"isolate":iso_,"V4_copy":ci,"V4_120bp":s120}
        # meta-analysis
        if s90 in ST.index:
            x=ST.loc[s90]; r.update(meta_prevalence_pct=round(x.prevalence_pct,2),meta_mean_rel_pct=round(x.mean_rel_pct,4),meta_max_rel_pct=round(x.max_rel_pct,3),meta_rank=int(rk.get(s90,np.nan)) if s90 in rk.index else None,
                meta_percentile=round(100*rk[s90]/len(present),2) if s90 in rk.index else None,meta_taxonomy=f"{str(tax.loc[s90].Hit_name).split(' 16S')[0]} ({tax.loc[s90].Per_ident}%)" if s90 in tax.index else "")
        else: r.update(meta_prevalence_pct=0.0,meta_mean_rel_pct=0.0,meta_note="V4 not present in meta-analysis")
        # AmphiBac isolates sharing this V4 (exact 90 bp contained)
        ab=[k for k,v in fa.items() if s90 in v or s90 in v.translate(rc)[::-1]]
        f_=m.set_index("FastaID").loc[ab,"Bd_function"].value_counts().to_dict() if ab else {}
        r.update(AmphiBac_isolates_same_V4=len(ab),AmphiBac_Bd_function=str(f_),AmphiBac_hosts="; ".join(m.set_index("FastaID").loc[ab,"Host.Species"].value_counts().index[:4]) if ab else "")
        for lab_,(t_,q_) in {"bullfrog":(bft,bfq),"salamander":(sat,saq)}.items():
            s=asv_stats(t_,q_,s120); r.update({f"{lab_}_match":s["match"],f"{lab_}_prevalence_pct":s["prev"],f"{lab_}_mean_rel_pct":s["mean"],f"{lab_}_max_rel_pct":s["max"]})
        out.append(r)
ISO=pd.DataFrame(out)
pd.set_option("display.width",300,"display.max_colwidth",60)
print(ISO.drop(columns=["V4_120bp"]).to_string())
