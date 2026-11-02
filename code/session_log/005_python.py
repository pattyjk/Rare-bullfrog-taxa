# session cell 5 | 2026-10-15 19:40 | env: python | language: python
# files written: []


for k in ["BF4","BF5"]:
    print("####",k); print(full[k][:1800])
f=re.sub(r"\s+"," ",full["BF2"])
print("####BF2", f[f.find("Distances"):f.find("Genome relatedness")+20])
