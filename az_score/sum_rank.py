import pandas as pd 

df = pd.read_csv("./rank.xls",sep="\t",index_col="CFD bin")
print(df.sum().sum())