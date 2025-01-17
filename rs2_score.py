from pathlib import Path
import pandas as pd
import numpy as np
from pysam import FastxFile
from Bio.Seq import Seq
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


# flank up down stream 10bp and calc Azimuth score
def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def reference_obtain(nc_no: str) -> str:
    # read ref
    ref = f"/mnt/ntc_data/wayne/Repositories/CRISPR/GCF_000001405.40/split_fa/{nc_no}.fa"
    ref_seq = ""
    for seq in FastxFile(ref):
        ref_seq = seq.sequence
    return ref_seq


def flanking_seq(ref_seq: str, start: int, end: int) -> list:
    left_terminal_pos = min(start, end)
    right_terminal_pos = max(start, end)
    
    left_flanking_seq = ref_seq[left_terminal_pos - 11 : left_terminal_pos - 1]
    right_flankig_seq = ref_seq[right_terminal_pos : right_terminal_pos + 10]
    
    if end > start:
        return [left_flanking_seq.upper(), right_flankig_seq.upper()]
    else:
        # reverse complement and switch left,right
        return [reverse_complement(right_flankig_seq).upper(), reverse_complement(left_flanking_seq).upper()]


def az_score(nc_no,seq):
    import azimuth.model_comparison
    # predictions = azimuth.model_comparison.predict(seq)
    # split numpy array
    batch_size = 2000
    scores = np.array([])
    for i in range(0,len(seq),batch_size):
        batch = seq[i:i+batch_size]
        score = azimuth.model_comparison.predict(seq=batch,model_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/score/Azimuth/azimuth/saved_models/V3_model_nopos.pickle")
        print(f"<{nc_no}> {i+batch_size+1}/{len(seq)} finished...")
        scores = np.append(scores,score)
    return scores


def rs2_score_calc(nc_no: str, gdb_df: pd.DataFrame) -> None:
    # concat full seq to calc rs2 score
    full_seq = gdb_df.apply(
    lambda row: row["up_stream"][-4:] + row[1] + row["down_stream"][:6],
    axis=1
)   
    full_seq = full_seq.to_numpy()
    score_li = az_score(nc_no,full_seq)
    gdb_df['az_score'] = score_li


def run(nc_no: str) -> None:
    # if result already exists,skiped
    if Path(f"/mnt/ntc_data/wayne/Repositories/CRISPR/az_score/{nc_no}.tsv").exists():
        print(f"/mnt/ntc_data/wayne/Repositories/CRISPR/az_score/{nc_no}.tsv is exists,skiped !!!")
        return None
    
    flank_res_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/flank_fill/{nc_no}.tsv"
    gdb_df = pd.DataFrame([])
    # if no flanking
    if not Path(flank_res_path).exists():
        # gdb 2 df 
        gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/cfd_filter/{nc_no}.tsv"
        type_li =  ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category", "float", "int32"]
        gdb_df = tsv2df(gdb_path,type_li)
        print(f"{nc_no} gdb load finished...")
        # flanking up down stream
        ref_seq = reference_obtain(nc_no)
        gdb_df[['up_stream', 'down_stream']] = gdb_df.apply(lambda x:pd.Series(flanking_seq(ref_seq,x[5],x[6])),axis=1)
        output_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/flank_fill/{nc_no}.tsv"
        gdb_df.to_csv(output_path,sep="\t",header=None,index=None)
        print(f"{nc_no} flank fill finished...")
    
    else:
        # if flank have exists
        gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/flank_fill/{nc_no}.tsv"
        type_li =  ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category", "float", "int32", "string", "string"]
        gdb_df = tsv2df(gdb_path,type_li)
        # new_col_dict =
        gdb_df.rename(columns={
            20:"up_stream",
            21:"down_stream"
        },inplace=True)
        print(f"/mnt/ntc_data/wayne/Repositories/CRISPR/flank_fill/{nc_no} load ...")

    # calc rs2 score
    rs2_score_calc(nc_no, gdb_df)
    print(f"{nc_no} rs2 score finished...")
    
    # save as tsv
    output_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/az_score/{nc_no}.tsv"
    gdb_df.to_csv(output_path,sep="\t",header=None,index=None)
    print(f"{nc_no} done,saved in /mnt/ntc_data/wayne/Repositories/CRISPR/az_score/{nc_no}.tsv")


def main() -> None:
    nc2chr_file = "/mnt/ntc_data/wayne/Repositories/CRISPR/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    async_in_iterable_structure(run,nc_li,24)
    # run(nc_li[-4])
    
if __name__ == "__main__":
    main()