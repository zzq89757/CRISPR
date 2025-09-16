import pandas as pd
from collections import defaultdict
from sys import path
path.append("/mnt_data/Wayne/Repositories/CRISPR/")
from generate_split_ori import async_in_iterable_structure


def re_num(nc_no: str) -> None:
    df_li = []
    insertion_count = 0
    # 数据库路径
    old_ko_path = f"/mnt_data/Wayne/Repositories/CRISPR/mysql/split_gtf/{nc_no}.tsv"
    fresh_ko_path = f"/mnt_data/Wayne/Repositories/CRISPR/mysql/re_num/{nc_no}.tsv"
    db_ai_path = f"/mnt_data/Wayne/Repositories/CRISPR/database_ai/snp_mark/{nc_no}.tsv"
    re_num_path = f"/mnt_data/Wayne/Repositories/CRISPR/database_ai/re_num/{nc_no}.tsv"
    # 文件读取为df
    old_ko_df = pd.read_csv(old_ko_path,sep="\t")
    fresh_ko_df = pd.read_csv(fresh_ko_path,sep="\t",header=None)
    db_ai_df = pd.read_csv(db_ai_path,sep="\t",header=None)
    # 读取 ai 库 group by gene id
    for gene_id, sub_df in db_ai_df.groupby(8, sort=False):
        old_ko_only_dict = defaultdict(list)
        fresh_ko_only_dict = defaultdict(list)
        old_fresh_common_dict = defaultdict(list)
        sub_df.reset_index(drop=True, inplace=True)
        sub_df['raw_no'] = sub_df[0].str.split("_").str[-1].astype(int)
        old_ko_sub_df = old_ko_df[old_ko_df['gene_id'] == gene_id]
        fresh_ko_sub_df = fresh_ko_df[fresh_ko_df[10] == gene_id]
        # 若gene id对应的旧ko及新ko库均为空 直接使用index + 1为编号
        if len(old_ko_sub_df) == len(fresh_ko_sub_df) == 0:
            sub_df['new_name'] = sub_df[7] + "[gRNA" + (sub_df.index + 1).astype(str) + "]"
            new_header = ['new_name', 'raw_no'] + list(range(1,19))
            sub_df = sub_df[new_header]
            df_li.append(sub_df)
            # print(sub_df)
            # exit()
            continue
        # 对比三个库的序列 寻找新旧ko库交集、新旧独有 并使用ai库进行遍历
        # 旧ko->新ko grna_number和新版命名
        # ai->新ko 原始编号判断

        fresh_ko_sub_df['no'] = (
            fresh_ko_sub_df[0]
            .str.extract(r"\[gRNA(\d+)\]")  # 提取数字部分
        ).astype(int)
        insertion_set = set(fresh_ko_sub_df['no'].to_list()) & set(old_ko_sub_df["grna_number"].to_list())
        fresh_ko_sub_df['is_insertion'] = fresh_ko_sub_df['no'].isin(insertion_set)
        old_ko_sub_df['is_insertion'] = old_ko_sub_df['grna_number'].isin(insertion_set)
        # 先用fresh 全部 和ai找交集 再找old only 和ai交集并编号 
    
        # 记录max number
        max_fresh_num = max(fresh_ko_sub_df['no']) if len(fresh_ko_sub_df) else 0
        max_old_num = max(old_ko_sub_df['grna_number']) if len(old_ko_sub_df) else 0
        max_num = max(max_fresh_num, max_old_num)
        # 寻找新ko和ai交集 byseq
        raw_no_new_seq_map = dict(zip(fresh_ko_sub_df[2],fresh_ko_sub_df[0]))
        sub_df['new_name'] = sub_df[1].map(raw_no_new_seq_map)
        # 找old only 和ai交集
        old_ko_only_sub_df = old_ko_sub_df[old_ko_sub_df['is_insertion']==False]
        old_ko_only_sub_df['guide'] = old_ko_only_sub_df['seq'].str[:20]
        old_seq2number = dict(zip(old_ko_only_sub_df['guide'],old_ko_only_sub_df['grna_number']))
        sub_df['new_name'] = sub_df.apply(lambda x: old_seq2number.get(x[1],None) if pd.isna(x['new_name']) else x['new_name'], axis=1)
        sub_df['new_name'] = sub_df.apply(lambda x:x[7] + "[gRNA" + str(max_num + x.name + 1) + "]" if pd.isna(x['new_name']) else x['new_name'],axis = 1)
        # 重设表头
        new_header = ['new_name', 'raw_no'] + list(range(1,19))
        sub_df = sub_df[new_header]
        df_li.append(sub_df)
        # print(sub_df)
        # print(old_ko_sub_df)
        # print(fresh_ko_sub_df)
        # print(old_ko_only_sub_df)
        # exit()
    res_df = pd.concat(df_li)
    res_df.to_csv(re_num_path,sep="\t",header=None,index=False)
    # print(res_df)

def main() -> None:
    nc2chr_file = "/mnt_data/Wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # re_num(nc_li[11])
    async_in_iterable_structure(re_num,nc_li,24)


if __name__ == "__main__":
    main()