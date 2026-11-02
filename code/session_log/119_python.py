# session cell 119 | 2026-10-18 11:29 | env: python | language: python
# files written: []


import pandas as pd, numpy as np
from scipy import stats as sst
def vid(fn): return host.artifacts(filename=fn,exact=True)["artifacts"][0]["latest_version_id"]
wb=host.artifact_path(vid("AmphiBac_meta_abundance.xlsx")); nm=host.artifact_path(vid("AmphiBac_near_matches.xlsx"))
iso=pd.read_excel(wb,sheet_name="AmphiBac_isolates",dtype={"Proportional.Growth.Bd":str})
so=pd.read_excel(wb,sheet_name="AmphiBac_sOTUs")
nd=pd.read_excel(nm,sheet_name="unmatched_V4_isolates")
cfc=pd.read_excel(nm,sheet_name="cultured_by_class")
print(iso.shape, so.shape, nd.shape); print(cfc.to_string())
meta=pd.read_csv(host.artifact_path("e9e6f252-b2f8-43f7-8080-ccdeecc5c8ca"),sep="\t",dtype=str,low_memory=False)
print(meta["Assay.Method"].value_counts().head(8).to_dict())
pg=pd.to_numeric(meta["Proportional.Growth.Bd"],errors="coerce"); print("numeric PG:",pg.notna().sum(), pg.describe().round(3).to_dict())
print(pd.crosstab(meta["Assay.Method"].fillna("NA"),pg.notna()).sort_values(True,ascending=False).head(8))
