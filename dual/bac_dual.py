from pathlib import Path
import time
import pandas as pd
from sys import path
from os import system
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure
pd.options.mode.copy_on_write = True


def extract_pair_grna_info(pair_list: list) -> list:
    # pair list 为二维数组 每个元素为成对df
    res_li = []
    for pair in pair_list:
        grna1 = pair[0]
        grna2 = pair[1]
        # 提取方向列计算切点偏移量
        grna1_offset = grna1[5].astype(str) + "0.5"
        grna2_offset = grna2[5].astype(str) + "0.5"
        # 记录公共信息Gene ID、gene name
        gene_id = grna1[10].values[0]
        gene_name = grna1[9].values[0]
        # 计算距离
        distance = int(abs((grna1[8] + grna1_offset.astype(float)).values[0] -
                       (grna2[8] + grna2_offset.astype(float)).values[0]))

        # 分别提取每条grna的id、序列、pam、flank、location、tran cov、strand orientation、score、snp
        info_idx = [0, 1, 21, 2, 3, 23, 13, 17, 5, 19, 22, 25]
        merge_df = pd.concat(
            [grna1[info_idx].iloc[0], grna2[info_idx].iloc[0]]) if grna1[1].values[0] < grna2[1].values[0] else pd.concat(
            [grna2[info_idx].iloc[0], grna1[info_idx].iloc[0]])
        merge_df[28] = gene_id
        merge_df[29] = gene_name
        merge_df[30] = distance
        merge_df = merge_df.reset_index(drop=True)
        merge_df = pd.DataFrame(merge_df).T
        merge_df = merge_df[[24, 25, 26] + list(range(24))]
        res_li.append(merge_df)
    return res_li


def add_pair_id() -> None:
    ...


def low_ranking(low_df: pd.DataFrame) -> pd.DataFrame:
    ...


def dual(raw_db: str) -> pd.DataFrame:
    # 读取filter20 数据库
    df = tsv2df(raw_db, [])
    all_res_li = []
    # 按照基因分组
    for gene, sub_df in df.groupby(9, sort=False):
        # print(gene,end="\n")
        grna_num = len(sub_df)
        # 若只有一条 无法成对 跳过
        if grna_num == 1:
            continue
        # 若只有两条 考虑方向后配对并输出信息 并且pari id为1
        if grna_num == 2:
            grna1_df = sub_df.iloc[[0]]
            grna2_df = sub_df.iloc[[1]]
            if grna1_df[5].values[0] == grna2_df[5].values[0]:
                continue
            all_res_li += extract_pair_grna_info([[grna1_df, grna2_df]])
            continue
        # 两条以上的情况 暂时有放回
        sub_df = sub_df.reset_index(drop=True)
        # 先按照方向分组
        df_pos = sub_df[sub_df[5] == '+']
        df_neg = sub_df[sub_df[5] == '-']
        low_pair_li = []
        filter_pair_li = []
        for _, grna1 in df_pos.iterrows():
            # 若过滤前组合数小于20 输出原始排列结果 否则筛选
            if len(df_pos) * len(df_neg) > 20:
                # 两个分数差值都在20以内
                cfd_score = grna1[19]
                cfd_low = cfd_score - 20
                cfd_high = cfd_score + 20
                rs2_score = grna1[22]
                rs2_low = rs2_score - 20
                rs2_high = rs2_score + 20
                filter_candidate_df = df_neg[(df_neg[19] <= cfd_high) & (df_neg[19] >= cfd_low) & (
                    df_neg[22] <= rs2_high) & (df_neg[22] >= rs2_low)]
                # 两条gRNA中点间隔在50-1wbp
                grna_start = grna1[6]
                grna_end = grna1[7]
                grna_middle = (grna_end + grna_start) / 2
                middle_start_left = grna_middle - 5000
                middle_end_left = grna_middle - 50
                middle_start_right = grna_middle + 50
                middle_end_right = grna_middle + 5000
                filter_candidate_df['grna_middle'] = (
                    filter_candidate_df[6] + filter_candidate_df[7]) / 2
                filter_candidate_df = filter_candidate_df[
                    (middle_start_left < filter_candidate_df['grna_middle']) & (filter_candidate_df['grna_middle'] < middle_end_left) |
                    (middle_start_right < filter_candidate_df['grna_middle']) & (
                        filter_candidate_df['grna_middle'] < middle_end_right)
                ]
                del filter_candidate_df['grna_middle']
            # 将正向grna和负向过滤候选加入filter pair li
            for __, grna2 in filter_candidate_df.iterrows():
                grna1_df = pd.DataFrame(grna1).T
                grna2_df = pd.DataFrame(grna2).T
                filter_pair_li.append([grna1_df, grna2_df])
            # 将正向grna和负向low候选加入low pair li
            low_df = df_neg[~df_neg[1].isin(filter_candidate_df[1])]
            for __, grna2 in low_df.iterrows():
                grna1_df = pd.DataFrame(grna1).T
                grna2_df = pd.DataFrame(grna2).T
                low_pair_li.append([grna1_df, grna2_df])
        # 对low进行划分等级并排序

        # 处理final pair li 过滤后结果不足20条则从low_sort补充至20条 filter distance从小到大
        filter_pair_li.sort(key=lambda x:abs((x[1][8].values[0] + (x[1][5].astype(str) + "0.5").astype(float).values[0])-(x[0][8].values[0])+ (x[0][5].astype(str) + "0.5").astype(float)).values[0])
        
        final_pair_li = filter_pair_li[:20] if len(
            filter_pair_li) >= 20 else filter_pair_li + low_pair_li[:20 - len(filter_pair_li)]

        all_res_li += extract_pair_grna_info(final_pair_li)
    # 8589908  8590507  9140859  9141958
    # 拼接结果 添加pair id
    all_pair_df = pd.concat(all_res_li)
    print(all_pair_df)

def run_dual(nc_no: str) -> None:
    t1 = time.time()
    raw_db = f"/mnt/ntc_data/wayne/Repositories/CRISPR/filter_20/{nc_no}.tsv"
    output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/dual/{nc_no}.tsv"
    if Path(output).exists():
        print(f"{nc_no} exists !!!")
        return
    new_df = dual(raw_db)
    print(f"dual time cost:{time.time() - t1}")
    return
    # 保存为tsv
    new_df.to_csv(output, sep="\t", header=None, index=None)


def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # async_in_iterable_structure(run_dual, nc_li, 24)
    run_dual(nc_li[-1])


if __name__ == "__main__":
    main()
