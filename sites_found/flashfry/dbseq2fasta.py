import pandas as pd

df = pd.read_csv("../../exon_filter/NC_000024.10.tsv",sep="\t",header=None)

output_fa = open("y_sg.fasta",'w')

seq = df[1].to_numpy() + df[2].to_numpy()
seq_li = [">"] * len(df[0]) + df[8].to_numpy() + ["_"] * len(df[0]) + df[0].astype(str).to_numpy() + ["\n"] * len(df[0]) + seq

print("\n".join(seq_li),file=output_fa)