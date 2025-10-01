from glob import glob
import pandas as pd



def remove_n0(project_dir: str, chr2nc_dict:dict) -> None:
    # 先根据raw id 去重 再根据序列找到重复行对应的索引存入array 再根据array中的索引 去除行
    gdb_path = f"{project_dir}/ag_mark/"
    file_li = glob(f"{gdb_path}/NC*tsv")
    all_db_li = []
    for file in file_li:
        gdb = pd.read_csv(file,sep="\t",header=None)
        all_db_li.append(gdb)
        # print(f"{file} append !!!")
    all_gdb_df = pd.concat(all_db_li)

    # 根据 raw_id 去重
    rm_overlap_df = all_gdb_df.drop_duplicates(subset=[0])
    # 根据序列找到序列相同、id不同的重复行的id
    duplicates_id_li = rm_overlap_df[rm_overlap_df.duplicated(subset=[1],keep=False)][0].to_list()
    # duplicates = rm_overlap_df[rm_overlap_df.duplicated(subset=[1],keep=False)]
    # duplicates.to_csv("test.tsv",header=None,index=False,sep="\t")
    # print(duplicates)
    # print(all_gdb_df)
    all_gdb_df = all_gdb_df[~all_gdb_df[0].isin(duplicates_id_li)]
    # 按照染色体拆分生成文件
    for chr, sub_df in all_gdb_df.groupby(3,sort=False):
        nc_no = chr2nc_dict[chr]
        # 保存为文件
        output_file = f"{project_dir}/ag_mark_n0_rm/{nc_no}.tsv"
        sub_df.to_csv(output_file,header=None,index=False,sep="\t")


def main() -> None:
    nc_li = pd.read_csv("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/nc_li", header=None, sep="\t")[0].to_list()
    chr_li = pd.read_csv("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/nc_li", header=None, sep="\t")[1].to_list()
    nc2chr_dict = dict(zip(nc_li, chr_li))
    chr2nc_dict = dict(zip(chr_li, nc_li))
    remove_n0("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus", chr2nc_dict)


if __name__ == "__main__":
    main()