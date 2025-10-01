import pandas as pd


def flashfry_score_file2tsv(score_file: str) -> pd.DataFrame:
    score_df = pd.read_csv(score_file,sep="\t",usecols=["contig","DoenchCFD_specificityscore","otCount"])
    return score_df


def append_cfd_score(project_dir: str, nc_no: str) -> None:
    # 读取flashfry score table
    score_file = f"{project_dir}/flashfry_out/score_out/{nc_no}.tsv"
    score_df = flashfry_score_file2tsv(score_file)
    score_df.columns = [0,1,2]
    score_df.sort_values(0,inplace=True)
    
    # 读取gdb
    gdb_file = f"{project_dir}/ag_mark_n0_rm/{nc_no}.tsv"
    gdf = pd.read_csv(gdb_file,sep="\t",header=None)
    # 外连接（outer join）
    merged_df = pd.merge(gdf, score_df, on=0, how='outer')
    # 去除Nan
    merged_df.dropna(inplace=True)
    # 将ot count 调整为int类型
    merged_df['2_y'] = merged_df['2_y'].astype(int)
    # 导出文件
    output_file = f"{project_dir}/cfd_score/{nc_no}.tsv"
    merged_df.to_csv(output_file,sep="\t",index=False,header=None)


def main() -> None:
    append_cfd_score("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus","NC_000087.8")


if __name__ == "__main__":
    main()
