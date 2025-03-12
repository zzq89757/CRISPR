from pathlib import Path
import time
import pandas as pd
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure
pd.options.mode.copy_on_write = True


def distance_cal(grna1_df: pd.Series, grna2_df: pd.Series) -> int:
    # 切点距离间隔 （含方向偏移量校正）
    distance = abs(int((grna1_df[8] + float(grna1_df[5] + "0.5")) - (grna2_df[8] + float(grna2_df[5] + "0.5"))))  
    return distance


def distance_rank(distance: int) -> int:
    # 将distance分为4个范围，30bp-5kb；5kb-10kb；10kb-50kb；>50kb 由于基因长度普遍在几kb 因此从小到大判断节省时间 30bp以内的优先级高于50kb以外
    if distance < 30:return 3
    if distance < 5000:return 0
    if distance < 10000:return 1
    if distance < 50000:return 2
    return 4


def sum_score(grna1_df: pd.Series, grna2_df: pd.Series) -> float:
    return grna1_df[19] + grna2_df[19] + grna1_df[22] + grna2_df[22]


def filter_mark(distance: int, grna1_df: pd.Series, grna2_df: pd.Series) -> int:
    # 距离间隔判断 间隔在30bp - 3kb之间
    if distance < 30 or distance > 3000:
        return 0
    # 分数差值判断
    cfd_diff = abs(grna1_df[19] - grna2_df[19])
    rs2_diff = abs(grna1_df[22] - grna2_df[22])
    if cfd_diff > 20 or rs2_diff > 20:
        return 0
    return 1


def transform_index_pair_li(idx_pair_li: list, gene_df: pd.DataFrame) -> pd.DataFrame:
    all_pair_df_li = []
    for idx_pair in idx_pair_li:
        idx1, idx2, raw_pair_id, distance, filter_flag = idx_pair
        grna1 = gene_df.loc[idx1]
        grna2 = gene_df.loc[idx2]
        gene_id = grna1[10]
        chr_no = grna1[4]
        gene_name = grna1[9]
        pair_id = f"{gene_name}[pair#{raw_pair_id}]"
        is_cd = "n" if grna1[5] == grna2[5] else "y"
        # 分别提取每条grna的name(0),rawid(1)、upstream(21)、seq(2)、pam(3)、downstream(23)、orientation_chr(5)、cutsite_chr(8)、loc_gene(13)、tran cov(17)、orientation_gene(12)、cfd score(19)、rs2 score(22)、snp(25)
        info_idx = [0, 1, 21, 2, 3, 23, 5, 8, 13, 17, 12, 19, 22, 25]
        # 拼接两行 添加pair_ID distance gene_ID chr gene_name
        merge_ser = pd.concat([pd.Series([pair_id, distance, gene_id, chr_no, gene_name, is_cd]), grna1[info_idx], grna2[info_idx]], ignore_index=True)
        # 转换为df并转置
        merge_df = pd.DataFrame(merge_ser).T
        all_pair_df_li.append(merge_df)
    # 合并所有子结果 并按照distance从小到大排列，如果distance一样，按gRNA Pair ID从小到大排列
    all_pair_df = pd.concat(all_pair_df_li, ignore_index=True)
    return all_pair_df


def dual(raw_db: str) -> pd.DataFrame:
    # 读取filter20 数据库
    df = tsv2df(raw_db, [])
    filter_df_li = []
    # 按照基因分组
    for gene, sub_df in df.groupby(9, sort=False):
        # print(gene,end="\t")
        grna_num = len(sub_df)
        # 若只有一条 无法成对 跳过
        if grna_num == 1:
            continue
        # 两条及以上的情况 暂时有放回
        sub_df = sub_df.reset_index(drop=True)
        # 直接穷举所有pair 按照min rawID 从小到大排序
        all_idx_pair_li = []
        index_len = len(sub_df.index)
        for idx1 in range(index_len):          
            grna_pos = sub_df.loc[idx1]
            pos_rawID = grna_pos[1]
            for idx2 in range(idx1 + 1,index_len):
                grna_neg = sub_df.loc[idx2]
                neg_rawID = grna_neg[1]
                # 计算distance
                distance = distance_cal(grna_pos, grna_neg)
                # 过滤器筛查
                filter_passed = filter_mark(distance, grna_pos, grna_neg)
                # 比较正负向rawID大小 并将正负索引值按照其rawID从小到大排序再加上min rawID以及distance和filter标记(默认为0 通过过滤器则为1)后存入数组
                tmp_li = [idx1, idx2, pos_rawID, distance, filter_passed] if pos_rawID < neg_rawID else [idx2, idx1, neg_rawID, distance, filter_passed]
                all_idx_pair_li.append(tmp_li)
        # 对all_idx_pair_li按照min rawID从小到大进行排序 索引号+1即为原始pairID
        all_idx_pair_li.sort(key=lambda x:x[2])
        # 将all_idx_pair_li中的minrawID改为pairID号(索引+1)
        for idx in range(len(all_idx_pair_li)):
            all_idx_pair_li[idx][2] = idx + 1
        # 将all_idx_pair_li的数据进行信息拼接并添加pair id
        res_df = transform_index_pair_li(all_idx_pair_li, sub_df)     
        # 按照过滤器分离all_idx_pair_li 并将filter li 按照距离排序
        filtered_idx_pair_li = [x for x in all_idx_pair_li if x[-1]]
        filtered_idx_pair_li.sort(key=lambda x:x[3])
        # 若filter已足够 则取前二十个
        final_idx_pair_li = filtered_idx_pair_li[:20]
        # 从unfilter 回补
        if len(filtered_idx_pair_li) < 20:
            unfiltered_idx_pair_li = [x for x in all_idx_pair_li if x[-1] == 0]
            # 对unfiltered_idx_pair_li 按照距离分级(50bp-10kb；10kb-50kb；>50kb)以及总分从高到低排序
            sorted_unfiltered_idx_pair_li = sorted(unfiltered_idx_pair_li, key=lambda x: (distance_rank(x[3]), -sum_score(sub_df.loc[x[0]], sub_df.loc[x[1]])))
            # 去unfiltered回补
            final_idx_pair_li += unfiltered_idx_pair_li[:20-len(filtered_idx_pair_li)]
        # final_idx_pair_li 先按distance从小到大排列，再按gRNA Pair ID从小到大排列
        sorted_final_idx_pair_li = sorted(final_idx_pair_li, key= lambda x: (x[3], x[2]))
        # 将sorted_final_idx_pair_li的数据进行信息拼接并添加pair id
        res_df = transform_index_pair_li(sorted_final_idx_pair_li, sub_df)
        filter_df_li.append(res_df)
    # 合并所有子结果
    filter_df = pd.concat(filter_df_li, ignore_index=True)
    return filter_df
        
    # # 8589908  8590507  9140859  9141958

def run_dual(nc_no: str) -> None:
    t1 = time.time()
    raw_db = f"/mnt/ntc_data/wayne/Repositories/CRISPR/filter_50/{nc_no}.tsv"
    output = f"/mnt/ntc_data/wayne/Repositories/CRISPR/dual_cd/{nc_no}.tsv"
    if Path(output).exists():
        print(f"{nc_no} exists !!!")
        return
    new_df = dual(raw_db)
    # 保存为tsv
    new_df.to_csv(output, sep="\t", header=None, index=None)
    print(f"<{nc_no}> dual time cost:{round(time.time() - t1, 2)}s, dual_20 line count:{len(new_df)}")


def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # async_in_iterable_structure(run_dual, nc_li, 24)
    run_dual(nc_li[-1])


if __name__ == "__main__":
    main()
