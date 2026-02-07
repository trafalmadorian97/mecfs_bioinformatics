from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.axis
from matplotlib.axes import Axes

from mecfs_bio.util.plotting.save_fig import write_plots_to_dir

"""
Purpose is generate plots that illustrate the concepts
underlying linkage disequilibrium score regression
"""

_ld_score_label = "ld_score"
_h2 = "h2"

def gen_fake_data_df(slope:float,noise_level: float, ld_score_values:np.ndarray, seed: int ) -> pd.DataFrame:
    gen =np.random.default_rng(seed=seed)
    h2_values = np.maximum(slope*ld_score_values+1+ gen.normal(0,noise_level,size=len(ld_score_values)),

                           0)
    return pd.DataFrame(
        {
    _ld_score_label:ld_score_values,
         _h2: h2_values,
        }

    )



def draw_schematic_ldsc_on_ax(ax: Axes, slope: float,seed:int, title:str|None=None):
    df = gen_fake_data_df(slope=slope,
                          noise_level=0.1,
                          ld_score_values=np.linspace(0.2,1,100),
                          seed=seed
                          )
    ax.scatter(df[_ld_score_label], df[_h2])
    ax.plot([0.23,0.97],[0.23*slope+1, slope*0.97+1] ,linewidth=3, color="black")
    ax.set_xlim(0.19,1.01)
    ax.set_ylim(0.5,2)
    ax.tick_params(labelbottom=False, labelleft=False)
    ax.set_xlabel("LD Score")
    ax.set_ylabel("Wald Test Statistic")
    if title is not None:
        ax.set_title(title)

def draw_basic_ldsc_schematic():
    fig, axes= plt.subplots(1, 2,figsize=(10,10))
    draw_schematic_ldsc_on_ax(axes[0], slope=0.2, title="Low $h^2$", seed=100)
    draw_schematic_ldsc_on_ax(axes[1], slope=0.8, title=   "High $h^2$", seed=101 )
    fig.show()
    figs = {"ldsc_schematic":fig}
    write_plots_to_dir(Path("data"), figs)
    print("yo")

if __name__ == '__main__':
    draw_basic_ldsc_schematic()
