from os import system

from sys import path
import time

import pandas as pd
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR")

from process_border_withid import async_in_iterable_structure


def run(nc_no: str) -> None:
    t1 = time.time()
    # nc_no = "NC_000024.10"
    

    nc_table = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    df = pd.read_csv(nc_table, sep="\t", header=None)
    nc2chr_dict = df.set_index(0)[1].to_dict()
    raw_db_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/split_out/sorted/no_head/with_id/spCas9_Homo_{nc2chr_dict[nc_no]}.tsv"
    new_db_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/database_ai/nag_remove/{nc_no}.tsv"
    system(f"awk '$3~/GG$/' {raw_db_path} > {new_db_path}")
    

def main() -> None:
    nc_table = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_li = pd.read_csv(nc_table, sep="\t", header=None)[0].to_list()
    async_in_iterable_structure(run,nc_li,24)
    # run(nc_li[-1])

if __name__ == "__main__":
    main()