import pandas as pd
raw_df = pd.read_csv("raw_0.sli",sep="\t",header=None)
ft_df = pd.read_csv("filter_0.sli",sep="\t",header=None)

gene_df = pd.read_csv("/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/extract/Gene_list.tsv",sep="\t",header=None)
gene_dict = dict(zip(gene_df[0],gene_df[5]))


# raw_df = pd.merge(raw_df,gene_df,on=0)[[0,5]]
# ft_df = pd.merge(ft_df,gene_df,on=0)[[0,5]]
raw_df[1] = raw_df[0].apply(lambda x:gene_dict[x])
ft_df[1] = ft_df[0].apply(lambda x:gene_dict[x])

raw_df.to_csv("raw_0.nli",sep="\t",header=None,index=None)
ft_df.to_csv("filter_0.nli",sep="\t",header=None,index=None)