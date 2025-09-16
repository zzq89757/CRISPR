import time
import pandas as pd


def annotation_gdb(gdb: pd.DataFrame, gene_pos_df: pd.DataFrame) -> pd.DataFrame:
    # 需要分别考虑切点和起止 三个位置 超出基因末端时 loc 采用(+/- len(over))
    gene_pos_len = len(gene_pos_df)
    current_gene_pos_idx = 0
    for idx, g_item in gdb.iterrows():
        # print(current_gene_pos_idx)

        # 跳过染色体前端未到达基因区的gRNA以及索引切换后位于基因起始前的gRNA
        if g_item[6] < gene_pos_df[1][current_gene_pos_idx]:
            continue
        
        if g_item[6] > gene_pos_df[2][current_gene_pos_idx]:
            current_gene_pos_idx += 1
            # 索引超出时退出循环
            if current_gene_pos_idx == gene_pos_len:
                break
        if (
            gene_pos_df[1][current_gene_pos_idx]
            <= g_item[6]
            <= gene_pos_df[2][current_gene_pos_idx]
        ):
            # 记录信息
            ...
            print(f"{g_item[6]} in {gene_pos_df.iloc[current_gene_pos_idx][1]}-{gene_pos_df.iloc[current_gene_pos_idx][2]} ")
            # print(f"{g_item[6]} in {gene_pos_df.iloc[current_gene_pos_idx][1]}-{gene_pos_df.iloc[current_gene_pos_idx][2]} ", flush=True, end="\r")
            
    
    return pd.DataFrame([])


def main() -> None:
    t = time.time()
    nc_no = "NC_000024.10"
    chr = 'chrY'
    
    nc_table = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    df = pd.read_csv(nc_table,sep="\t",header=None)
    chr2nc_dict = df.set_index(1)[0].to_dict()
    type_li = ["string", "int32", "int32", "category", "string", "category"]
    type_dict = dict(enumerate(type_li))
    # 读取基因位置信息文件
    gene_pos_df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_gtf/{chr2nc_dict[chr]}/Gene_list_cut_insertion.tsv", sep="\t", header=None, low_memory=False, dtype=type_dict)
    # 读取gRNA db 文件 
    type_li = ['string', 'string', 'category', 'category', 'int32', 'int32', 'int32']
    type_dict= dict(enumerate(type_li))
    # gdb_df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/spCas9_Homo_{chr}.tsv", header=None, sep="\t")
    gdb_df = pd.read_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/y_10w.tsv", header=None, sep="\t")
    # gdb 注释
    annotated_df = annotation_gdb(gdb_df, gene_pos_df)
    # 保存为tsv文件
    # annotated_df.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/{chr}_gene.tsv", index=False, header=None, sep="\t")
    print(f"time cost:{time.time() - t}")
    

if __name__ == "__main__":
    main()
    