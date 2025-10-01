import pandas as pd
from pysam import FastxFile
from Bio.Seq import Seq
from .filter_intron import gdb2df

# flank up down stream 10bp and calc Azimuth score
def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def reference_obtain(project_dir: str, nc_no: str) -> str:
    # read ref
    ref = f"{project_dir}/GCF/fa/{nc_no}.fa"
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



def flank_fill(project_dir: str, nc_no: str) -> None:
    # gdb 2 df 
    gdb_path = f"{project_dir}/cfd_score/{nc_no}.tsv"
    gdb_df = gdb2df(gdb_path)
    # flanking up down stream
    ref_seq = reference_obtain(project_dir,nc_no)
    gdb_df[['up_stream', 'down_stream']] = gdb_df.apply(lambda x:pd.Series(flanking_seq(ref_seq,x[5],x[6])),axis=1)
    # filter flanking contains N
    gdb_df = gdb_df[~(gdb_df['up_stream'].str.contains('N') | gdb_df['down_stream'].str.contains('N'))]
    flank_res_path = f"{project_dir}/flank_fill/{nc_no}.tsv"
    gdb_df.to_csv(flank_res_path,sep="\t",header=None,index=None)
    


def main() -> None:
    nc2chr_file = "{project_dir}/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # async_in_iterable_structure(run,nc_li,24)
    flank_fill("NC_000017.11")
    
if __name__ == "__main__":
    main()