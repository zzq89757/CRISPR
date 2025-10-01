import pandas as pd
import numpy as np
from .filter_intron import gdb2df


def az_score(nc_no,seq):
    import azimuth.model_comparison
    # predictions = azimuth.model_comparison.predict(seq)
    # split numpy array
    batch_size = 2000
    scores = np.array([])
    for i in range(0,len(seq),batch_size):
        batch = seq[i:i+batch_size]
        score = azimuth.model_comparison.predict(seq=batch,model_file = "{project_dir}/score/Azimuth/azimuth/saved_models/V3_model_nopos.pickle")
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


def rs2_score(project_dir: str, nc_no: str) -> None:
    # flank db read
    gdb_path = f"{project_dir}/flank_fill/{nc_no}.tsv"
    gdb_df = gdb2df(gdb_path)
    # new_col_dict
    gdb_df.rename(columns={
        20:"up_stream",
        21:"down_stream"
    },inplace=True)
    # calc rs2 score
    rs2_score_calc(nc_no, gdb_df)
    # save as tsv
    output_path = f"{project_dir}/az_score/{nc_no}.tsv"
    gdb_df.to_csv(output_path,sep="\t",header=None,index=None)


def main() -> None:
    nc2chr_file = "{project_dir}/nc2chr.tsv"
    nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    nc_li = nc_df[0].tolist()
    # async_in_iterable_structure(run,nc_li,24)
    rs2_score("NC_000017.11")
    
if __name__ == "__main__":
    main()