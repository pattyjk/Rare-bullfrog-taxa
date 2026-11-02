# session cell 103 | 2026-10-17 18:12 | env: python | language: python
# files written: []


z=(PSm.cultured_sOTUs_3mm==0); print("samples:",len(PSm),"with zero matches:",int(z.sum()), "zero exact:",int((PSm.cultured_sOTUs_exact==0).sum()))
print(PSm.cultured_sOTUs_3mm.describe().round(1).to_dict())
print("spearman matches vs %reads:",sst.spearmanr(PSm.cultured_sOTUs_3mm,PSm.pct_reads_cultured_3mm))
print("spearman richness vs matches:",sst.spearmanr(PSm.richness,PSm.cultured_sOTUs_3mm).statistic)
print(PSm[z][["richness","total_reads","host_species","country"]].head(10).to_string())
