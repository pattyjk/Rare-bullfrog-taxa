# session cell 88 | 2026-10-17 17:51 | env: python | language: python
# files written: []


pres=set(a.index)
st_=iso.sOTUs_all.apply(lambda s:any(x in pres for x in s.split(";")))
print("Bd-tested matched isolates:",len(iso),"| with sOTU detected in amphibian data:",int(st_.sum()), iso[st_].Bd_function.value_counts().to_dict(), "| only absent:", int((~st_).sum()), iso[~st_].Bd_function.value_counts().to_dict())
print("frac inhibitory detected vs absent:", round((iso[st_].Bd_function=="Inhibitory").mean(),3), round((iso[~st_].Bd_function=="Inhibitory").mean(),3))
print(sst.chi2_contingency(pd.crosstab(st_, iso.Bd_function=="Inhibitory"))[:2])
