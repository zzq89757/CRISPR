
import pandas as pd

    
def gdb2df(gdb_path: str, type_li: list) -> pd.DataFrame:
    # type_li = ["int32", "string", "category", "category", "category", "int32", "int32", "int32", "string", "int32", "category", "category", "string", "int32"]
    type_dict = dict(enumerate(type_li))

    gdb_df = pd.read_csv(
        gdb_path,
        sep="\t",
        header=None,
        dtype=type_dict,
    )
    return gdb_df