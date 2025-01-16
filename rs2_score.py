from pathlib import Path
import pandas as pd
from pysam import FastxFile
from Bio.Seq import Seq
from sys import path
path.append("/mnt/ntc_data/wayne/Repositories/CRISPR/")
from utils.read_tsv import tsv2df
from generate_split_ori import async_in_iterable_structure


# 补齐上下游10bp后计算Azimuth score
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


def az_score(seq):
    import azimuth.model_comparison
    predictions = azimuth.model_comparison.predict(seq)
    return predictions


def rs2_score_calc(gdb_df: pd.DataFrame) -> None:
    # concat full seq to calc rs2 score
    full_seq = gdb_df.apply(
    lambda row: row["up_stream"][-4:] + row[1] + row["down_stream"][:6],
    axis=1
)   
    full_seq = full_seq.to_numpy()
    score_li = az_score(full_seq)
    gdb_df['az_score'] = score_li


def run(nc_no: str) -> None:
    # gdb 2 df 
    gdb_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/cfd_filter/{nc_no}.tsv"
    type_li =  ["string", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32", "string", "category", "string", "category", "float", "int32"]
    gdb_df = tsv2df(gdb_path,type_li)
    
    # flanking up down stream
    ref_seq = reference_obtain(nc_no)
    gdb_df[['up_stream', 'down_stream']] = gdb_df.apply(lambda x:pd.Series(flanking_seq(ref_seq,x[5],x[6])),axis=1)
    
    # calc rs2 score
    rs2_score_calc(gdb_df)
    
    # save as tsv
    output_path = f"/mnt/ntc_data/wayne/Repositories/CRISPR/az_score/{nc_no}.tsv"
    gdb_df.to_csv(output_path,sep="\t",header=None,index=None)


def main() -> None:
    run("NC_000024.10")
    
    
if __name__ == "__main__":
    main()