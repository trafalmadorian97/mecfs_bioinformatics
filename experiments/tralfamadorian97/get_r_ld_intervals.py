import numpy as np
import tarfile
from pathlib import Path
import pyreadr
import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects import numpy2ri
from rpy2.robjects.conversion import localconverter
def go():
    df = pd.read_csv("r_ld_file_list.txt",sep=r"\s+",
                     names=[
                         f"col_{i}" for i in range(6)
                     ])
    split = df["col_5"].str.split(".",expand=True)
    bounds = split[3].str.split("_",expand=True)
    bounds = bounds.astype(int)
    pos =   173846152
    import pdb; pdb.set_trace()
    filtered = bounds.loc[(bounds.iloc[:,0]<= pos)&(bounds.iloc[:,1]>= pos)]
    import pdb; pdb.set_trace()
    print("yo")

def print_members():
    src_path = Path("1.tar.gz")
    with tarfile.open(src_path, "r:gz") as tar_object:
        for member in tar_object.getmembers():
            print(member.name)


def retrieve():
   src_path = Path("1.tar.gz")
   with tarfile.open(src_path, "r:gz") as tar_object:
       tar_object.extract("1/ukb_b38_0.1_chr1.R_snp.173128768_175120632.Rvar",
                          "ukb_b38_0.1_chr1.R_snp.173128768_175120632.Rvar")

       tar_object.extract("1/ukb_b38_0.1_chr1.R_snp.173128768_175120632.RDS",
                          "ukb_b38_0.1_chr1.R_snp.173128768_175120632.RDS")

def read_rds_matrix(rds_path: str | Path):
    res = pyreadr.read_r(str(rds_path))  # returns an OrderedDict-like
    if len(res.keys()) != 1:
        raise ValueError(f"Expected 1 object in RDS, found: {list(res.keys())}")
    obj = next(iter(res.values()))

    # obj is often a pandas DataFrame (if it was an R matrix/data.frame)
    return obj

def read_rds_dense_as_numpy(rds_path: str | Path) -> np.ndarray:
    readRDS = ro.r["readRDS"]
    x = readRDS(str(rds_path))

    # If it's a base R matrix, as.matrix() works
    x_mat = ro.r["as.matrix"](x)

    with localconverter(ro.default_converter + numpy2ri.converter):
        arr = np.array(x_mat)
    return arr

def get_ld_matrix():
    pth = Path("ukb_b38_0.1_chr1.R_snp.173128768_175120632.RDS/1/ukb_b38_0.1_chr1.R_snp.173128768_175120632.RDS")
    obj = read_rds_dense_as_numpy(pth)
    import pdb; pdb.set_trace()
    print("hello")


if __name__ == "__main__":
    get_ld_matrix()
    # print_members()
    # retrieve()


#            0          1
#0   173128768  175120632
#75  173128768  175120632