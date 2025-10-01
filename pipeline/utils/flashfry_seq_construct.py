import pandas as pd


def construct_seq(project_dir: str, nc_no: str) -> None:
    gdb_path = f"{project_dir}/ag_mark_n0_rm/{nc_no}.tsv"
    output_seq_file = f"{project_dir}/flashfry_out/fa/{nc_no}.fa"
    gdf = pd.read_csv(gdb_path,header=None,sep="\t",usecols=[0,1,2]).drop_duplicates()
    tmp_str = "".join([f">{id}\n{seq}{pam}\n" for id, seq, pam in zip(gdf[0],gdf[1],gdf[2])])
    output_handle = open(output_seq_file,'w')
    output_handle.write(tmp_str)


def main() -> None:
    construct_seq("/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus","NC_000086.8")


if __name__ == "__main__":
    main()