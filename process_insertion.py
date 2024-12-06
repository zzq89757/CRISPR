from instertion_find import process_gene_pos_insertion, finnaly_check
import pandas as pd
from collections import Counter

pd.options.mode.copy_on_write = True


# 拆分函数
def split_intervals(df: pd.DataFrame, cuts: list) -> pd.DataFrame:
    cuts = sorted(cuts)  # 切点排序
    result = []

    for _, row in df.iterrows():
        gene, start, end, strand = row[0], row[1], row[2], row[3]
        segment_start = start
        cut_flag = 0
        for cut in cuts:
            # 可切区间
            if cut > segment_start and cut <= end:
                cut_flag = 1
                result.append([gene, segment_start + 1, cut - 1, strand])
                segment_start = cut
        # 不可切区间
        if not cut_flag:
            result.append([gene, start + 1, end - 1, strand])

    # print(pd.DataFrame(result))
    return pd.DataFrame(result)


# 合并基因子区间
def merge_intervals(df: pd.DataFrame) -> pd.DataFrame:
    merged_li = []
    for common_start, sub_df in df.groupby(1):
        sub_df = sub_df.reset_index(drop=True)
        merged_item = [
            ",".join(sub_df[0].to_list()),
            common_start,
            sub_df[2][0],
            ",".join(sub_df[3].to_list()),
        ]
        merged_li.append(merged_item)
    # print(pd.DataFrame(merged_li))
    return pd.DataFrame(merged_li)


# 回补切点
def fill_cut_site(
    regions_df: pd.DataFrame, merged_sub_df: pd.DataFrame, cut_site_li: list
) -> pd.DataFrame:
    # 获取切点覆盖基因列表
    cover_gene_li = []
    filled_cut_site_li = [pd.DataFrame([])] * len(merged_sub_df)
    filled_region_idx_li = []
    for cut_site in cut_site_li:
        sub_df = regions_df[(regions_df[1] <= cut_site) & (cut_site <= regions_df[2])]
        cover_gene_li.append(sub_df[0].to_list())
    # print(cover_gene_li)
    # print(merged_sub_df)
    # 回补切点
    for i, cut_site in enumerate(cut_site_li):
        sub_df = merged_sub_df[
            (merged_sub_df[1] == cut_site + 1) | (cut_site == merged_sub_df[2] + 1)
        ]
        sub_df = sub_df[sub_df[0] == ",".join(cover_gene_li[i])]
        # 获取待插入切点区间的索引号
        idx = int(sub_df.index[0])
        filled_region_idx_li.append(idx)

        tmp_sub_df = (
            filled_cut_site_li[idx]
            if not filled_cut_site_li[idx].empty
            else sub_df.loc[idx]
        )
        # tmp_sub_df = filled_cut_site_li[idx]
        # print(tmp_sub_df)

        # 填充区间
        if tmp_sub_df[1] == cut_site + 1 or tmp_sub_df[1] == cut_site - 1:
            tmp_sub_df[1] = cut_site
        if tmp_sub_df[2] == cut_site - 1 or tmp_sub_df[2] == cut_site + 1:
            tmp_sub_df[2] = cut_site
        filled_cut_site_li[idx] = tmp_sub_df
        # print(cut_site)
        # print(sub_df)
    # print(cut_site_li)
    # print(filled_cut_site_li)
    # print(merged_sub_df)
    new_merged_li = [
        merged_sub_df.iloc[i] if len(v) == 0 else v
        for i, v in enumerate(filled_cut_site_li)
    ]
    # print(new_merged_li)
    # print(pd.DataFrame(new_merged_li))
    return pd.DataFrame(new_merged_li)


def insertion_merge(regions_li: list) -> pd.DataFrame:
    regions_df = pd.DataFrame(regions_li)
    regions_df = regions_df.reset_index(drop=True)
    # print(regions_df)
    # 根据所有区间起止拆分子区间后将各个基因放入

    # 检查子区间内基因数目 根据每个子区间交集的起始（切点）设置位置偏移量

    # 切点包括所有的区间起始终止 (排除min start 和max end)
    min_start = regions_df[1][0]
    max_end = max(regions_df[2])
    # print(min_start,end="\t")
    # print(max_end)
    cut_site_li = sorted(pd.concat([regions_df[1], regions_df[2]]).drop_duplicates())
    # 根据切点拆分子区间
    cutted_df = split_intervals(regions_df, cut_site_li)
    # 合并子区间相同的item
    merged_df = merge_intervals(cutted_df)
    # 看切点落在哪几个基因区间里 后续进行回补
    filled_df = fill_cut_site(regions_df, merged_df, cut_site_li)
    return filled_df
    # print(cut_site_li)

    # exit()


def insertion_detective(gene_pos_table) -> pd.DataFrame:
    type_li = ["string", "int32", "int32", "category"]
    type_dict = dict(enumerate(type_li))
    res = []
    # 示例输入数据
    df = pd.read_csv(
        gene_pos_table, sep="\t", header=None, low_memory=False, dtype=type_dict
    )
    tmp_res = []
    common_idx_li = []
    for i in range(len(df) - 1):
        # 跳过交集区
        if i in common_idx_li:
            continue
        # 无重叠 直接append
        if df.iloc[i + 1][1] > df.iloc[i][2]:
            res.append(df.iloc[i])
        else:
            # 记录MAX end 并将后续区间同maxend比较以判断是否为连续交集区
            max_end = df.iloc[i][2]
            for j in range(i, len(df) - 1):
                if df.iloc[j][2] > max_end:
                    max_end = df.iloc[j][2]
                tmp_res.append(df.iloc[j])
                common_idx_li.append(j)
                # print(f"max_end is {max_end}")
                # 下一区间的起始大于maxend 退出循环
                if df.iloc[j + 1][1] > max_end:
                    # print(j)
                    break
                # 判断最后一行是否可并入连续交集区
                if j == len(df) - 2:
                    if df.iloc[j][1] < max_end:
                        tmp_res.append(df.iloc[j + 1])
                        common_idx_li.append(j + 1)
            # tmp res 用于处理公共区的结果
            # print(common_idx_li)
            for i in insertion_merge(tmp_res).iterrows():
                # print(i[1])
                res.append(i[1])

            tmp_res = []
    result_df = pd.DataFrame(res)

    return result_df


def check_all():
    from pathlib import Path

    for pt in Path("split_gtf/").glob("*"):
        gtf_file = pt / "Gene_list.tsv"
        finnaly_check(gtf_file)


def main() -> None:
    # merged_df = process_gene_pos_insertion("split_gtf/NC_000001.11/Gene_list.tsv")
    # check_all()
    # merged_df = insertion_detective("/mnt/ntc_data/wayne/Repositories/CRISPR/test.tsv")
    merged_df = insertion_detective("/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/NC_000001.11/Gene_list.tsv")
    merged_df.to_csv("./test_g.tsv", index=False, header=None, sep="\t")


if __name__ == "__main__":
    main()
