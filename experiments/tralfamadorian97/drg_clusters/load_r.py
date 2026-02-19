import pyreadr
import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

from rpy2.robjects import numpy2ri, pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import (
    importr,
)
import rpy2.robjects as ro

def go():
    result = pyreadr.read_r("human_meta_final_cluster.Rdata")
    import pdb; pdb.set_trace()
    print(result.keys())

def go25():
    result = pyreadr.read_r("Hs_LCM_Seurat_export_slim.RData")
    import pdb; pdb.set_trace()
    print(result.keys())

def go2():
    result = pyreadr.read_r("Hs_LCM_final.RData")
    # import pdb;
    # pdb.set_trace()
    print(result.keys())

def go3():

    # 1. Load the R instance and libraries
    seurat = importr('Seurat')
    base = importr('base')

    # 2. Load the .RData file
    # This loads the object into the embedded R environment
    conv = ro.default_converter + pandas2ri.converter + numpy2ri.converter
    robjects.r['load']("Hs_LCM_final.RData")
    counts = robjects.r["HS.counts.raw"]
    meta = robjects.r["Meta.fused"]
    meta_2 = robjects.r["Meta.1066.wo.glia"]
    with localconverter(conv):
        counts_py = ro.conversion.get_conversion().rpy2py(counts)
        meta_py=ro.conversion.get_conversion().rpy2py(meta)
        meta_2_py = ro.conversion.get_conversion().rpy2py(meta_2)
    import pdb; pdb.set_trace()
    print("yo")

if __name__ == '__main__':
    go3()