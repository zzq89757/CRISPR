import pandas as pd
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df


nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
nc_li = nc_df[0].tolist()

all_res = pd.DataFrame([])
# for nc_no in nc_li:
for nc_no in nc_li:
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/rs2_score/{nc_no}.tsv"
    # gdb_path = "/mnt/ntc_data/wayne/Repositories/CRISPR/ag_mark/NC_000024.10.tsv"
    # read gdb
    header_type_li = ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category"]
    gdb_df = tsv2df(gdb_path, [])
    # 定义区间
    bins_cfd = [0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
    bins_rule = [0.01, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00]
    # 添加区间列
    gdb_df['CFD bin'] = pd.cut(gdb_df[13], bins=bins_cfd, right=True)
    gdb_df['RS bin'] = pd.cut(gdb_df[17], bins=bins_rule, right=True)
    
    # 生成交叉表
    result = pd.crosstab(gdb_df['CFD bin'], gdb_df['RS bin'])
    # 格式化显示
    result.index = result.index.astype(str)
    result.columns = result.columns.astype(str)

    # 输出结果
    # print(result)
    if nc_no == "NC_000001.11":
        all_res = result
    else:
        all_res += result
all_res.to_csv("/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/two_score_rank-Andy.xls",sep="\t")