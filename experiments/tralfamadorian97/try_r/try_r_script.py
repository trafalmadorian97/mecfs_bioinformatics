import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


# R package names

# R vector of strings
from rpy2.robjects.vectors import StrVector
import rpy2.robjects as robjects

from rpy2.robjects import pandas2ri
import rpy2.robjects.packages as rpackages
from rpy2.robjects.packages import importr


packnames = ('ggplot2', 'hexbin', 'BiocManager')

def install_r_packages():
    # import R's utility package
    utils = rpackages.importr('utils')

    # select a mirror for R packages
    utils.chooseCRANmirror(ind=1)  # select the first mirror in the list
    # Selectively install what needs to be install.
    # We are fancy, just because we can.
    names_to_install = [x for x in packnames if not rpackages.isinstalled(x)]
    if len(names_to_install) > 0:
        utils.install_packages(StrVector(names_to_install))

def try_celldex():
    manager=  rpackages.importr('BiocManager')

    manager.install("celldex")
    se = importr('SummarizedExperiment')
    celldex = importr('celldex')
    base = importr("base")
    celldex_reference =  celldex.fetchReference("immgen", version="2024-02-26")
    r_col_data = se.colData(celldex_reference)
    r_df = base.as_data_frame(r_col_data)
    with (ro.default_converter + pandas2ri.converter).context():
        df = ro.conversion.get_conversion().rpy2py(r_df)
    # sample_metadata_df = pandas2ri.rpy2_to_pandas_dataframe(r_col_data)
    rlist =list(base.rownames(r_col_data))
    import pdb; pdb.set_trace()
    import pdb; pdb.set_trace()
    print("yo")
    # rpackages.importr("bioconductor-celldex")





def go():
    print("yo")


if __name__ == "__main__":
    # install_r_packages()
    try_celldex()
    go()