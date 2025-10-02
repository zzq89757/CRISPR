import pandas as pd
import numpy as np

def gdb2df(gdb_path: str) -> pd.DataFrame:
    # type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32"]

    # type_dict = dict(enumerate(type_li))

    gdb_df = pd.read_csv(
        gdb_path,
        sep="\t",
        header=None,
        # dtype=type_dict,
        # low_memory=False,
    )
    # print(gdb_df)
    return gdb_df


def az_score(nc_no,seq):
    import azimuth.model_comparison
    # predictions = azimuth.model_comparison.predict(seq)
    # split numpy array
    batch_size = 2000
    scores = np.array([])
    for i in range(0,len(seq),batch_size):
        batch = seq[i:i+batch_size]
        score = azimuth.model_comparison.predict(seq=batch,model_file = "/mnt_data/Wayne/Software/miniconda3/envs/azimuth3/lib/python3.8/site-packages/Biomatters_Azimuth-0.2.6-py3.8.egg/azimuth/saved_models/V3_model_nopos.pickle")
        print(f"<{nc_no}> {i+batch_size}/{len(seq)} finished...")
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
    output_path = f"{project_dir}/rs2_score/{nc_no}.tsv"
    gdb_df.to_csv(output_path,sep="\t",header=None,index=None)


def main() -> None:
    from sys import argv
    # nc2chr_file = "{project_dir}/nc2chr.tsv"
    # nc_df = pd.read_csv(nc2chr_file, sep="\t", header=None)
    # nc_li = nc_df[0].tolist()
    # async_in_iterable_structure(run,nc_li,24)
    rs2_score(argv[1],argv[2])

    
if __name__ == "__main__":
    main()