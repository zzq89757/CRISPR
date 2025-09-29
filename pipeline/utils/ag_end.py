import pandas as pd
from .filter_intron import gdb2df


def mark_ag(gdb: pd.DataFrame) -> None:
    gdb[17] = gdb[1].apply(lambda x: x[-2:] if str(x).endswith(('AG', 'GG')) else 'no')


def ag_mark(project_dir: str, nc_no: str) -> None:
    gdb_path = f"{project_dir}/tran_count/{nc_no}.tsv"
    gdb = gdb2df(gdb_path)
    # 注释gdb
    mark_ag(gdb)
    gdb.to_csv(f"{project_dir}/ag_mark/{nc_no}.tsv",header=None,sep="\t",index=False)
    

def main() -> None:
    nc2chr_file = "/mnt_data/Wayne/Repositories/CRISPR/pipeline/mus/nc_li"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    ag_mark(nc_li[-1])
    
    
if __name__ == "__main__":
    main()
    