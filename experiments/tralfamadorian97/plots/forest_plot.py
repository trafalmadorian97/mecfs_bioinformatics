from pathlib import Path

import matplotlib.pyplot as plt
import forestplot as fp

from mecfs_bio.util.plotting.save_fig import write_plots_to_dir


def go():

    df = fp.load_data("sleep")  # companion example data
    fig,ax=plt.subplots(1,1,figsize=(16, 8))
    fp.forestplot(df,  # the dataframe with results data
                  estimate="r",  # col containing estimated effect size
                  ll="ll", hl="hl",  # columns containing conf. int. lower and higher limits
                  varlabel="label",  # column containing variable label
                  ylabel="Confidence interval",  # y-label title
                  xlabel="Pearson correlation",  # x-label title
                  ax=ax
                  )
    figs = {
        "forest":fig
    }
    write_plots_to_dir(path=Path("data"),plots=figs)

if __name__ == '__main__':
    go()

