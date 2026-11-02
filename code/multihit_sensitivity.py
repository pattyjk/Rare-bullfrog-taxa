#!/usr/bin/env python3
"""Multi-hit sensitivity of the AmphiBac abundance-class x Bd-function results.

Question: AmphiBac isolates and meta-analysis sOTUs match many-to-many (an isolate can be tied
between several nearest sOTUs; an sOTU can carry isolates with conflicting Bd labels). Does the
way these multi-hits are resolved change the dominant / intermediate / rare results?

Run from the folder containing inputs/:  python multihit_sensitivity.py
Inputs : inputs/amphib_only_otu_table.zip (meta-analysis 90-bp V4 sOTU x sample counts)
         inputs/AmphibBac_2025.1-tested.fasta, inputs/AmphibBac_Metadata_2025.1-tested.tsv
Outputs: amphibac_isolate_sOTU_candidates.csv.gz, meta_sOTU_stats.csv.gz,
         AmphiBac_multihit_sensitivity.png, AmphiBac_multihit_sensitivity.xlsx
Requires: python>=3.10, numpy, pandas, scipy, matplotlib, openpyxl. Peak memory ~4.5 GB; runtime ~2 min.
Assembled verbatim from the analysis session cells (diagnostic prints kept), with input paths
made local and the plotting-style helpers replaced by the simple stand-ins below.
"""
import matplotlib; matplotlib.use("Agg")
from collections import Counter
import matplotlib.pyplot as _plt
def apply_figure_style(sizes=(8,7,6)):
    _plt.rcParams.update({"font.size":sizes[0],"axes.titlesize":sizes[0],"axes.labelsize":sizes[0],
        "xtick.labelsize":sizes[2],"ytick.labelsize":sizes[2],"legend.fontsize":sizes[1],
        "axes.spines.top":False,"axes.spines.right":False})
def panel_letter(ax,letter):
    ax.text(-0.02,1.08,letter,transform=ax.transAxes,fontweight="bold",fontsize=10,ha="right",va="bottom")
# ---- session cell 177 ----
import pandas as pd, numpy as np, re, zipfile, io, os, time
t0=time.time()
os.makedirs('work',exist_ok=True)
z=zipfile.ZipFile('inputs/amphib_only_otu_table.zip'); nm=[n for n in z.namelist() if n.endswith('.txt')][0]
z.extract(nm,'work'); txt='work/'+nm
with open(txt) as f: hdr=f.readline().rstrip('\n').split('\t')
nf=len(hdr)+1
ids=[];rows=[]
with open(txt) as f:
    f.readline()
    for line in f:
        p=line.rstrip('\n').split('\t')
        if len(p)!=nf or len(p[0])!=90: continue
        ids.append(p[0]); rows.append(np.array(p[1:],dtype=np.float32))
X=np.vstack(rows); del rows
tot=X.sum(0); keep=tot>=1000
R=X[:,keep]/tot[keep]*100
mean_rel=R.mean(1); prev=(R>0).mean(1)*100
st=pd.DataFrame({'sOTU':ids,'mean_rel_pct':mean_rel,'prevalence_pct':prev})
st=st[st.mean_rel_pct>0].reset_index(drop=True)
print(keep.sum(), len(ids), len(st), round(time.time()-t0))

# ---- session cell 178 ----
t0=time.time()
def read_fasta(p):
    d={};k=None;b=[]
    for l in open(p):
        l=l.strip()
        if l.startswith('>'):
            if k: d[k]=''.join(b)
            k=l[1:].split()[0];b=[]
        else: b.append(l.upper())
    if k: d[k]=''.join(b)
    return d
fa=read_fasta('inputs/AmphibBac_2025.1-tested.fasta')
meta=pd.read_csv('inputs/AmphibBac_Metadata_2025.1-tested.tsv',sep='\t',low_memory=False)
comp=str.maketrans('ACGTRYKMBDHVN','TGCAYRMKVHDBN')
rc=lambda s:s.translate(comp)[::-1]
P=re.compile('GTG[CT]CAGC[AC]GCCGCGGTAA')
def win(s,n=90):
    for q in (s,rc(s)):
        m=P.search(q)
        if m and m.end()+n<=len(q): return q[m.end():m.end()+n]
    return None
S=set(st.sOTU); ab=dict(zip(st.sOTU,st.mean_rel_pct))
M=np.frombuffer(''.join(st.sOTU).encode(),dtype=np.uint8).reshape(-1,90)
rec=[]
for k,s in fa.items():
    ex=set()
    for q in (s,rc(s)):
        for i in range(len(q)-89):
            w=q[i:i+90]
            if w in S: ex.add(w)
    w=win(s); dmin=None; near=[]
    if not ex and w is not None and re.fullmatch('[ACGTN]+',w):
        d=(M!=np.frombuffer(w.encode(),dtype=np.uint8)).sum(1); dmin=int(d.min())
        if dmin<=3: near=list(st.sOTU.values[d==dmin])
    rec.append((k,w is not None,';'.join(sorted(ex)),len(ex),dmin,';'.join(near),len(near)))
I=pd.DataFrame(rec,columns=['FastaID','has_V4','exact_sOTUs','n_exact','nearest_mm','near_sOTUs','n_near'])
print(len(I), I.has_V4.sum(), ((I.n_exact>0)|(I.n_near>0)).sum(), ((I.n_exact>0)|(I.n_near>0))[I.has_V4].sum(), round(time.time()-t0))

# ---- session cell 180 ----
I['has_V4b']=I.has_V4|(I.n_exact>0)
det=(I.n_exact>0)|(I.n_near>0)
print(I.has_V4b.sum(), det.sum(), (~I.has_V4b).sum())
bd=meta.set_index('FastaID')
print([c for c in meta.columns if 'Bd' in c or 'Ref' in c or 'Growth' in c])

# ---- session cell 182 ----
amb=I[(I.n_exact==0)&I.has_V4].copy()
amb['w']=[win(fa[k]) for k in amb.FastaID]
amb=amb[~amb.w.str.fullmatch('[ACGTN]+')]
d3=[int((M!=np.frombuffer(w.encode(),dtype=np.uint8)).sum(1).min()) for w in amb.w]
print(len(amb), sum(x<=3 for x in d3))

# ---- session cell 183 ----
import scipy.stats as ss
for idx in amb.index:
    w=amb.at[idx,'w']; d=(M!=np.frombuffer(w.encode(),dtype=np.uint8)).sum(1); dm=int(d.min())
    I.at[idx,'nearest_mm']=dm
    if dm<=3: near=list(st.sOTU.values[d==dm]); I.at[idx,'near_sOTUs']=';'.join(near); I.at[idx,'n_near']=len(near)
def norm_bd(x):
    x=str(x).strip()
    if x in ("Inhibitory","inhibitory","Inhibitor"): return "Inhibitory"
    if x in ("NonInhibitory","NoEffect","No Effect","No effect","Non-inhibitor","ns"): return "Non-inhibitory"
    if x in ("Facilitating","enhancing"): return "Facilitating"
    return None
I=I.merge(meta[['FastaID','Bd.Inhibition','Proportional.Growth.Bd','Reference']],on='FastaID',how='left')
I['Bd_function']=I['Bd.Inhibition'].map(norm_bd)
I['inhib_score']=(1-pd.to_numeric(I['Proportional.Growth.Bd'],errors='coerce'))*100
# candidate set = all exact sOTUs, else all tied nearest
I['cand']=np.where(I.n_exact>0,I.exact_sOTUs,I.near_sOTUs)
I['n_cand']=np.where(I.n_exact>0,I.n_exact,I.n_near)
print('detected',(I.n_cand>0).sum())
def most(c,f): 
    L=c.split(';'); return f(L,key=lambda s:ab[s])
# original rule: exact -> all; near tie -> most abundant
I['set_orig']=[c if ne>0 else (most(c,max) if c else '') for c,ne in zip(I.cand,I.n_exact)]
CL=['Dominant','Intermediate','Rare']
cls=lambda a:'Dominant' if a>=1 else ('Intermediate' if a>=0.1 else 'Rare')
st['cls']=st.mean_rel_pct.map(cls); clsd=dict(zip(st.sOTU,st.cls))
def taxa(col,df=None,rule='majority'):
    df=I if df is None else df
    b=df[df.Bd_function.notna()&(df[col]!='')]
    E=b.assign(sOTU=b[col].str.split(';')).explode('sOTU')
    g=E.groupby('sOTU')
    T=pd.DataFrame({'n':g.FastaID.nunique(),'fi':g.Bd_function.apply(lambda x:(x=='Inhibitory').mean()),
                    'ff':g.Bd_function.apply(lambda x:(x=='Facilitating').mean()),'med':g.inhib_score.median()})
    T['cls']=T.index.map(clsd); return T,E
T,E=taxa('set_orig')
T['inh']=T.fi>=0.5
print(T.cls.value_counts().reindex(CL).tolist(), (100*T.groupby('cls').inh.mean()).reindex(CL).round(0).tolist())
print('isolates',E.FastaID.nunique())

# ---- session cell 184 ----
out=I[['FastaID','Reference','Bd_function','inhib_score','has_V4b','n_exact','exact_sOTUs','nearest_mm','n_near','near_sOTUs','n_cand','set_orig']].rename(columns={'has_V4b':'has_V4'})
out.to_csv('amphibac_isolate_sOTU_candidates.csv.gz',index=False)
st.to_csv('meta_sOTU_stats.csv.gz',index=False)
b=I[I.Bd_function.notna()&(I.n_cand>0)]
print('Bd-tested matched isolates',len(b))
print('exact multi',(b.n_exact>1).sum(),' exact single',(b.n_exact==1).sum(),' near single',((b.n_exact==0)&(b.n_near==1)).sum(),' near tied',((b.n_exact==0)&(b.n_near>1)).sum())
nt=b[(b.n_exact==0)&(b.n_near>1)]
print('near-tie size median/max',nt.n_near.median(),nt.n_near.max())
# class change most-abundant vs rarest
ca=[clsd[most(c,max)] for c in nt.cand]; cr=[clsd[most(c,min)] for c in nt.cand]
print(pd.crosstab(pd.Series(ca,name='most_abund'),pd.Series(cr,name='rarest')))
em=b[b.n_exact>1]; print('exact-multi class mixes', Counter(tuple(sorted(set(clsd[s] for s in c.split(';')))) for c in em.cand).most_common(5))
# sOTU label conflict
Eo=E; g=Eo.groupby('sOTU').Bd_function
nl=g.nunique(); n=g.size()
print('sOTUs with >=2 isolates',(n>=2).sum(),' conflicting',((n>=2)&(nl>1)).sum(), round(100*((n>=2)&(nl>1)).sum()/(n>=2).sum(),1))
print('conflict by class', (pd.DataFrame({'n':n,'nl':nl,'cls':T.cls})[n>=2].assign(c=lambda d:d.nl>1).groupby('cls').c.agg(['size','mean'])).reindex(CL).round(2).to_dict())

# ---- session cell 185 ----
rng=np.random.default_rng(1)
def summarise(T,lab):
    t=T.copy(); t['inh']=lab
    k=pd.crosstab(t.cls,t.inh).reindex(CL).fillna(0)
    p=ss.chi2_contingency(k)[1] if k.shape[1]==2 else np.nan
    pct=(100*t.groupby('cls').inh.mean()).reindex(CL)
    rho=ss.spearmanr(t.index.map(ab),t.inh.astype(float))[0]
    return dict(n_D=int((t.cls=='Dominant').sum()),n_I=int((t.cls=='Intermediate').sum()),n_R=int((t.cls=='Rare').sum()),
                inh_D=pct['Dominant'],inh_I=pct['Intermediate'],inh_R=pct['Rare'],chi2_p=p,spearman_abund_vs_inh=rho)
res=[]
def add(name,group,T,lab): d=summarise(T,lab); d.update(scenario=name,group=group); res.append(d)
# --- tie rules (majority label)
I['set_rare']=[c if ne>0 else (most(c,min) if c else '') for c,ne in zip(I.cand,I.n_exact)]
I['set_all']=I.cand
I['set_unamb']=np.where(I.n_cand==1,I.cand,'')
I['set_exact']=np.where(I.n_exact>0,I.cand,'')
for col,name in [('set_orig','Tied nearest -> most abundant sOTU (used so far)'),('set_rare','Tied nearest -> rarest sOTU'),
                 ('set_all','Tied nearest -> all tied sOTUs'),('set_unamb','Drop isolates with tied nearest sOTUs'),('set_exact','Exact matches only (0 mismatches)')]:
    Tt,_=taxa(col); add(name,'A. isolate-to-sOTU tie rule (majority label)',Tt,Tt.fi>=0.5)
# random tie draws
rr=[]
for it in range(500):
    I['set_rand']=[c if ne>0 else (rng.choice(c.split(';')) if c else '') for c,ne in zip(I.cand,I.n_exact)]
    Tt,_=taxa('set_rand'); rr.append(summarise(Tt,Tt.fi>=0.5))
rr=pd.DataFrame(rr)
d={k:rr[k].median() for k in rr}; d.update(scenario='Tied nearest -> random sOTU (median of 500 draws)',group='A. isolate-to-sOTU tie rule (majority label)'); res.append(d)
print(rr[['inh_D','inh_I','inh_R','chi2_p']].quantile([.025,.5,.975]).round(3))

# ---- session cell 186 ----
To,Eo=taxa('set_orig')
add('Majority of isolates inhibitory (used so far)','B. sOTU label rule (original tie rule)',To,To.fi>=0.5)
add('Any isolate inhibitory (AmphiBac subset-mode logic)','B. sOTU label rule (original tie rule)',To,To.fi>0)
add('All isolates inhibitory','B. sOTU label rule (original tie rule)',To,To.fi==1)
un=Eo.groupby('sOTU').Bd_function.nunique()==1
Tu=To[un.reindex(To.index).values]
add('Unanimous sOTUs only (conflicting sOTUs dropped)','B. sOTU label rule (original tie rule)',Tu,Tu.fi==1)
Tm=To[(Eo.groupby('sOTU').size()>=2).reindex(To.index).values]
add('Majority, sOTUs with >=2 isolates only','B. sOTU label rule (original tie rule)',Tm,Tm.fi>=0.5)
# top-hit analogue: one random isolate per sOTU
rr2=[]
grp=Eo.groupby('sOTU').Bd_function.apply(list)
for it in range(500):
    lab=pd.Series({s:(rng.choice(v)=='Inhibitory') for s,v in grp.items()}).reindex(To.index)
    rr2.append(summarise(To,lab))
rr2=pd.DataFrame(rr2)
d={k:rr2[k].median() for k in rr2}; d.update(scenario='One random isolate per sOTU (top-hit analogue; median of 500 draws)',group='B. sOTU label rule (original tie rule)'); res.append(d)
print(rr2[['inh_D','inh_I','inh_R','chi2_p']].quantile([.025,.5,.975]).round(3))
R_=pd.DataFrame(res)[['group','scenario','n_D','n_I','n_R','inh_D','inh_I','inh_R','chi2_p','spearman_abund_vs_inh']]
pd.set_option('display.width',250); print(R_.round(3).to_string(index=False))
# continuous score, each tie rule
sc=[]
for col,name in [('set_orig','most abundant'),('set_rare','rarest'),('set_all','all tied'),('set_unamb','drop tied'),('set_exact','exact only')]:
    Tt,_=taxa(col); Tt=Tt[Tt.med.notna()]
    sc.append(dict(tie_rule=name,median_D=Tt[Tt.cls=='Dominant'].med.median(),median_I=Tt[Tt.cls=='Intermediate'].med.median(),median_R=Tt[Tt.cls=='Rare'].med.median(),
                   kruskal_p=ss.kruskal(*[Tt[Tt.cls==c].med for c in CL])[1],spearman_abund_vs_score=ss.spearmanr(Tt.index.map(ab),Tt.med)[0],spearman_p=ss.spearmanr(Tt.index.map(ab),Tt.med)[1]))
SC=pd.DataFrame(sc); print(SC.round(3).to_string(index=False))

# ---- session cell 187 ----
print('isolates per sOTU (median, mean) by class:', To.groupby('cls').n.agg(['median','mean']).reindex(CL).round(1).to_dict())
iso=[];cul=[]
tot_cls=st.cls.value_counts().reindex(CL)
for col,name in [('set_orig','most abundant (used so far)'),('set_rare','rarest'),('set_all','all tied'),('set_unamb','drop tied'),('set_exact','exact only')]:
    b=I[I.Bd_function.notna()&(I[col]!='')].copy()
    if col=='set_all':  # isolate takes the class of its most abundant AND rarest? use fractional: each tied sOTU weight 1/k
        E2=b.assign(sOTU=b[col].str.split(';')).explode('sOTU'); E2['w']=1/E2.groupby('FastaID').sOTU.transform('size')
    else:
        E2=b.assign(sOTU=b[col]); E2['w']=1.0
    E2['cls']=E2.sOTU.map(clsd); E2['inh']=(E2.Bd_function=='Inhibitory')
    g=E2.groupby('cls').apply(lambda d:pd.Series({'n':d.w.sum(),'pct':100*(d.w*d.inh).sum()/d.w.sum()})).reindex(CL)
    iso.append(dict(tie_rule=name,n_D=g.n['Dominant'],n_I=g.n['Intermediate'],n_R=g.n['Rare'],inh_D=g.pct['Dominant'],inh_I=g.pct['Intermediate'],inh_R=g.pct['Rare']))
    a=I[I[col]!='']; cs=set(';'.join(a[col]).split(';'))
    c=st[st.sOTU.isin(cs)].cls.value_counts().reindex(CL).fillna(0)
    cul.append(dict(tie_rule=name,cultured_D=f"{int(c['Dominant'])}/{tot_cls['Dominant']}",cultured_I=f"{int(c['Intermediate'])}/{tot_cls['Intermediate']}",cultured_R=f"{int(c['Rare'])}/{tot_cls['Rare']}"))
ISO=pd.DataFrame(iso); CUL=pd.DataFrame(cul)
print(ISO.round(1).to_string(index=False)); print(CUL.to_string(index=False))

# ---- session cell 188 ----
import matplotlib as mpl, matplotlib.pyplot as plt
apply_figure_style(sizes=(8,7,6))
colc={'Dominant':'#08306B','Intermediate':'#4292C6','Rare':'#A6A6A6'}
short={'Tied nearest -> most abundant sOTU (used so far)':'Most abundant tied sOTU (used so far)',
 'Tied nearest -> rarest sOTU':'Rarest tied sOTU','Tied nearest -> all tied sOTUs':'All tied sOTUs',
 'Drop isolates with tied nearest sOTUs':'Drop tied isolates','Exact matches only (0 mismatches)':'Exact matches only',
 'Tied nearest -> random sOTU (median of 500 draws)':'Random tied sOTU (500 draws)',
 'Majority of isolates inhibitory (used so far)':'Majority inhibitory (used so far)',
 'Any isolate inhibitory (AmphiBac subset-mode logic)':'Any isolate inhibitory',
 'All isolates inhibitory':'All isolates inhibitory',
 'One random isolate per sOTU (top-hit analogue; median of 500 draws)':'One random isolate (top-hit; 500 draws)',
 'Unanimous sOTUs only (conflicting sOTUs dropped)':'Unanimous sOTUs only',
 'Majority, sOTUs with >=2 isolates only':'Majority, sOTUs with \u22652 isolates'}
fig,axs=plt.subplots(2,1,figsize=(6.6,4.6),sharex=True,gridspec_kw={'height_ratios':[6,6],'hspace':0.35})
for ax,grp,rand,title in [(axs[0],'A. isolate-to-sOTU tie rule (majority label)',rr,'Isolates tied between sOTUs: how they are assigned'),
                          (axs[1],'B. sOTU label rule (original tie rule)',rr2,'sOTUs with conflicting isolate labels: how the sOTU is called')]:
    d=R_[R_.group==grp].reset_index(drop=True)
    yy=np.arange(len(d))[::-1]
    for y,(_,r) in zip(yy,d.iterrows()):
        ax.plot([min(r.inh_D,r.inh_I,r.inh_R),max(r.inh_D,r.inh_I,r.inh_R)],[y,y],color='#DDDDDD',lw=1,zorder=1)
        israndom='500 draws' in r.scenario
        for c,k in [('Rare','inh_R'),('Intermediate','inh_I'),('Dominant','inh_D')]:
            if israndom:
                lo,hi=rand[k].quantile([.025,.975]); ax.plot([lo,hi],[y,y],color=colc[c],lw=2.2,alpha=0.45,solid_capstyle='butt',zorder=2)
            ax.scatter(r[k],y,s=26,color=colc[c],zorder=3,edgecolor='white',lw=0.4)
        p=r.chi2_p; ax.text(101,y,(f"p = {p:.2g}" if p>=0.001 else "p < 0.001")+(" (median)" if israndom else ""),va='center',fontsize=6,color='#333333')
    ax.set_yticks(yy); ax.set_yticklabels([short[s] for s in d.scenario])
    for t in ax.get_yticklabels():
        if 'used so far' in t.get_text(): t.set_fontweight('bold')
    ax.set_ylim(-0.6,len(d)-0.4); ax.set_title(title,loc='left')
axs[1].set_xlim(0,100); axs[1].set_xlabel('sOTUs called Bd-inhibitory (%)')
from matplotlib.lines import Line2D
h=[Line2D([],[],marker='o',ls='',color=colc[c],label=l,markersize=5) for c,l in [('Dominant','Dominant \u22651%'),('Intermediate','Intermediate 0.1\u20131%'),('Rare','Rare <0.1%')]]
axs[0].legend(handles=h,frameon=False,loc='lower left',bbox_to_anchor=(0.0,1.12),ncol=3,fontsize=6,handletextpad=0.2,columnspacing=0.8,borderaxespad=0)
axs[0].text(1.0,1.16,'\u03c7\u00b2 test, class \u00d7 inhibitory',transform=axs[0].transAxes,ha='right',fontsize=6,color='#333333')
for a_,l_ in zip(axs,'ab'): panel_letter(a_,l_)
fig.savefig('AmphiBac_multihit_sensitivity.png',dpi=300,bbox_inches='tight')
r_=fig.canvas.get_renderer()
texts=[(t,t.get_window_extent(r_)) for t in fig.findobj(mpl.text.Text) if t.get_text().strip() and t.get_visible()]
print([(a.get_text()[:25],b.get_text()[:25]) for i,(a,ba) in enumerate(texts) for b,bb in texts[i+1:] if ba.overlaps(bb)])


# ---- session cell 189 ----
b=I[I.Bd_function.notna()&(I.n_cand>0)]
prev_tab=pd.DataFrame([
 ('Bd-tested isolates matched to a meta-analysis sOTU (<=3 mm)',len(b)),
 ('  exact containment of exactly one sOTU',int((b.n_exact==1).sum())),
 ('  exact containment of >1 sOTU',int((b.n_exact>1).sum())),
 ('  no exact match; single nearest sOTU (1-3 mm)',int(((b.n_exact==0)&(b.n_near==1)).sum())),
 ('  no exact match; nearest distance tied between >1 sOTU',int(((b.n_exact==0)&(b.n_near>1)).sum())),
 ('  tie size, median (max)',f"{nt.n_near.median():.0f} ({nt.n_near.max()})"),
 ('  tied isolates whose class changes most-abundant -> rarest',229),
 ('sOTUs with >=2 Bd-tested isolates',int((n>=2).sum())),
 ('  with conflicting Bd labels',f"{((n>=2)&(nl>1)).sum()} ({100*((n>=2)&(nl>1)).sum()/(n>=2).sum():.1f}%)"),
 ('  conflicting: dominant / intermediate / rare','13/13 (100%) / 45/49 (92%) / 188/265 (71%)'),
 ('Bd-tested isolates per sOTU, median: dominant / intermediate / rare','19 / 15 / 2')],columns=['quantity','value'])
# verify conflict counts by class
cc=pd.DataFrame({'n':n,'nl':nl,'cls':T.cls})[n>=2]; print(cc.groupby('cls').apply(lambda d:f"{(d.nl>1).sum()}/{len(d)}"))
rand_ci=pd.concat({'random tied sOTU':rr[['inh_D','inh_I','inh_R','chi2_p']].quantile([.025,.5,.975]),
                   'one random isolate per sOTU':rr2[['inh_D','inh_I','inh_R','chi2_p']].quantile([.025,.5,.975])})
notes=pd.DataFrame({'note':[
 'Question: does the many-to-many matching between AmphiBac isolates and meta-analysis sOTUs (the issue examined in the AmphiBac robustness thread) change the rank/abundance-class results of this project?',
 'Data: AmphiBac 2025.1-tested (7,873 isolates); meta-analysis 90-bp V4 sOTUs with reads in 2,411 samples with >=1,000 reads (48,431 sOTUs). Abundance class from mean relative abundance across all 2,411 samples: dominant >=1%, intermediate 0.1-1%, rare <0.1%.',
 'Matching (unchanged): exact containment of a 90-bp sOTU anywhere in the isolate 16S (both strands); otherwise nearest sOTU by Hamming distance over the 90 bp after 515F, accepted if <=3 mismatches (IUPAC codes count as mismatches). This reproduces the earlier numbers exactly (6,610 detected isolates; 15/55/472 Bd-tested taxa; 67/49/39% inhibitory; 5,972 isolates).',
 'Part A varies how isolates tied between several nearest sOTUs are assigned (label rule fixed at majority). Part B varies how an sOTU with conflicting isolate labels is called (tie rule fixed at most-abundant). Random scenarios: 500 draws; medians in the table, 2.5-97.5% ranges in sheet random_draw_ranges.',
 'chi2_p: chi-square test of abundance class x inhibitory (taxon level, the conservative unit). spearman_abund_vs_inh: Spearman correlation between sOTU mean relative abundance and the inhibitory call (0/1).',
 'Any/All label rules depend on the number of isolates per sOTU (dominant median 19, rare median 2), so they are confounded with sampling depth: more isolates make "any inhibitory" more likely and "all inhibitory" less likely.',
 'Isolate-level: for "all tied" each tied sOTU gets weight 1/k. Cultured share: sOTUs with >=1 matched AmphiBac isolate (Bd-tested or not).',
]})


# ---- session cell 190 ----
prev_tab.loc[prev_tab.quantity.str.contains('dominant / intermediate / rare') & prev_tab.quantity.str.contains('conflicting'),'value']='13/13 (100%) / 45/49 (92%) / 189/265 (71%)'
with pd.ExcelWriter('AmphiBac_multihit_sensitivity.xlsx') as w:
    for df,sh,ix in [(notes,'notes',False),(prev_tab,'multihit_prevalence',False),(R_,'taxon_pct_inhibitory',False),(rand_ci,'random_draw_ranges',True),
                     (SC,'taxon_score_by_tie_rule',False),(ISO,'isolate_level_by_tie_rule',False),(CUL,'cultured_share_by_tie_rule',False)]:
        df.to_excel(w,sheet_name=sh,index=ix)
print(prev_tab.to_string(index=False))