import pandas as pd
import time


def process_gene_pos_insertion(gene_pos_table: str) -> pd.DataFrame:
    type_li = ['string', 'int32', 'int32','category']
    type_dict= dict(enumerate(type_li))
    gene_df = pd.read_csv(
        gene_pos_table, sep="\t", header=None, low_memory=False, dtype=type_dict
    )
    # 转换为元组
    # gene_info_tuple = list(
    #     gene_df.itertuples(index=False)
    # )  # 设置 index=False 不包括索引

    # 合并公共区
    for i in range(1, len(gene_df)):
        # 检查是否有重叠
        overlap_start = 0
        overlap_end = 0
        # 若为同一区域
        # if gene_df[1][i] == gene_df[1][i - 1] and 
        if gene_df[1][i] < gene_df[2][i - 1]:  # 存在交集
            # 记录交集区域起止
            overlap_start = gene_df[1][i]
            overlap_end = gene_df[2][i - 1]
            # 在表末尾插入新的交集行
            gene_df.loc[len(gene_df)] = [
                f"{gene_df[0][i-1]},{gene_df[0][i]}",
                overlap_start,
                overlap_end,
                f"{gene_df[3][i-1]},{gene_df[3][i]}",
            ]
            # 修改原始行 去除重叠区
            gene_df[2][i - 1] = overlap_start - 1
            gene_df[1][i] = overlap_end + 1
    # 重新排序
    # gene_df[4] = gene_df[4].replace(";", "", regex=True)
    gene_df = gene_df.sort_values(by=[1]).reset_index(drop=True)
    print(gene_df)
    return gene_df


def process_sg(gdb_file: str, gene_pos_df: pd.DataFrame) -> None:
    current_idx = 0
    
    type_li = ['string', 'string', 'category', 'category', 'int32', 'int32', 'int32']
    type_dict= dict(enumerate(type_li))
    # 打开gRNA数据库文件
    gdb_df = pd.read_csv(gdb_file, sep="\t", header=None, low_memory=False, dtype=type_dict)

    for i in range(len(gdb_df)):
        # 索引超出时退出循环
        if current_idx == len(gene_pos_df):
            break

        # 跳过染色体前端未到达基因区的gRNA以及索引切换后位于基因起始前的gRNA
        if gdb_df[6][i] < gene_pos_df[1][current_idx]:
            continue
        # 索引记录已经超出的区域
        for j in range(current_idx, len(gene_pos_df)):
            if gdb_df[6][i] > gene_pos_df[2][current_idx]:
                current_idx = j
            else:
                if (
                    gene_pos_df[1][current_idx]
                    <= gdb_df[6][i]
                    <= gene_pos_df[2][current_idx]
                ):
                    # 记录信息
                    ...
                    # print(f"{gdb_df.iloc[i][6]} in {gene_pos_df.iloc[current_idx][1]}-{gene_pos_df.iloc[current_idx][2]} ")
                else:
                    ...
                    # print(f"{gdb_df.iloc[i][6]} not in {gene_pos_df.iloc[current_idx][1]}-{gene_pos_df.iloc[current_idx][2]} ")
                    # mihari
                    #- 安阿 以伊 宇 衣江 *於
                    #k 加 幾 久 计 已
                    #s 佐散 之 *寸须 世 曾
                    #t 太多 知千 川 天 止
                    #n 奈 仁 奴 *祢 乃
                    #h *波八 比 不 *部 保
                    #*m  末万 美三 武 女 毛
                    #y 也     由     与
                    #*r 良 利 留流 礼 吕
                    #*w 和           逺乎
                    #*[n] 无尔
                break
        if gdb_df[6][i] > gene_pos_df[2][current_idx]:
            continue
        # print(f"current_idx is {current_idx}")
        # 切点在基因区域内 标注基因信息

def finnaly_check(gene_pos_table):
    type_li = ['string', 'int32', 'int32','category']
    type_dict= dict(enumerate(type_li))
    gene_df = pd.read_csv(
        gene_pos_table, sep="\t", header=None, low_memory=False, dtype=type_dict
    )
    final_start = gene_df.iloc[len(gene_df) -1][1]
    for i in gene_df[2][:-1]:
        if final_start < i :
            print(gene_pos_table)
            print(i)


def main() -> None:
    gene_pos_table = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/NC_000001.11/Gene_list.tsv"
    gdb_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/test_gene/sg_sample.tsv"
    gdb_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/y_100w.tsv"
    # gdb_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/spCas9_Homo_chrY_nohead.tsv"
    merged_gene_df = process_gene_pos_insertion(gene_pos_table)
    t = time.time()
    process_sg(gdb_file, merged_gene_df)
    print(f"time cost:{time.time() - t}")


if __name__ == "__main__":
    main()
