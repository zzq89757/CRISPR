import pandas as pd

df1 = pd.read_csv("/mnt/ntc_data/wayne/Repositories/CRISPR/cds_mark_re/NC_000024.10.tsv",header=None, usecols=[15],sep="\t")[15]
df2 = pd.read_csv("/mnt/ntc_data/wayne/Repositories/CRISPR/cds_mark_re/NC_000024.10_s.tsv",header=None, usecols=[15],sep="\t")[15]

idx = 0
for i1, i2 in zip(df1, df2):
    if i1 != i2:
        print(idx)
    idx += 1