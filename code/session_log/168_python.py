# session cell 168 | 2026-10-28 18:35 | env: python | language: python
# files written: []

import pandas as pd
m=pd.read_csv(host.artifact_path('fb82e995-1af5-43e1-9f42-be8c34e85602')); print(m.columns.tolist()); print(m.iloc[:, :8].to_string())
a=pd.read_csv(host.artifact_path('107585dd-e41f-4250-9007-7b68bb9b313f')); print(len(a)); print(a.iloc[:,0].value_counts().to_string())