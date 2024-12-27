from curses import noecho
import pandas as pd



def gdb2df(gdb_path: str) -> pd.DataFrame:
    type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "string"]

    type_dict = dict(enumerate(type_li))

    gdb_df = pd.read_csv(
        gdb_path,
        sep="\t",
        header=None,
        dtype=type_dict,
        # low_memory=False,
    )
    # print(gdb_df)
    return gdb_df


def mark_ag(gdb: pd.DataFrame) -> None:
    gdb[17] = gdb[1].apply(lambda x: 'yes' if str(x).endswith(('AG', 'GG')) else 'no')



def run_ag_mark(nc_no: str) -> None:
    gdb_path = "/mnt/ntc_data/wayne/Repositories/CRISPR/tran_count/NC_000024.10.tsv"
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/tran_count/{nc_no}.tsv"
    
    gdb = gdb2df(gdb_path)
    
    mark_ag(gdb)
    # gdb.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/ag_mark/NC_000024.10.tsv",header=None,sep="\t")
    gdb.to_csv(f"/mnt/ntc_data/wayne/Repositories/CRISPR/ag_mark/{nc_no}.tsv",header=None,sep="\t")
    

def main() -> None:
    from generate_split_ori import async_in_iterable_structure
    nc2chr_file = "nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run_ag_mark, nc_li, 24)
    
    
if __name__ == "__main__":
    main()