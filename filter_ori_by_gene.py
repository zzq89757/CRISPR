import pandas as pd


def annotation_gdb(gdb: pd.DataFrame, gene_pos_df: pd.DataFrame) -> pd.DataFrame:
    
    return pd.DataFrame([])


def main() -> None:
    # 读取gRNA db 文件 
    gdb_df = pd.read_csv()
    # 读取基因位置信息文件
    gene_pos_df = pd.read_csv()
    # gdb 注释
    annotated_df = annotation_gdb(gdb_df, gene_pos_df)
    # 保存为tsv文件
    annotated_df.to_csv()
    

if __name__ == "__main__":
    main()
    