import matplotlib.image as mpimg
import numpy as np
import matplotlib.pyplot as plt
import zepid
from zepid.graphics import EffectMeasurePlot


def go():
    labs = ["ACA(Isq=41.37% Tausq=0.146 pvalue=0.039 )",
            "ICA0(Isq=25.75% Tausq=0.092 pvalue=0.16 )",
            "ICA1(Isq=60.34% Tausq=0.121 pvalue=0.00 )",
            "ICAb(Isq=25.94% Tausq=0.083 pvalue=0.16 )",
            "ICAw(Isq=74.22% Tausq=0.465 pvalue=0.00 )"]
    measure = [2.09, 2.24, 1.79, 2.71, 1.97]
    lower = [1.49, 1.63, 1.33, 2.00, 1.25]
    upper = [2.92, 3.07, 2.42, 3.66, 3.11]
    p = EffectMeasurePlot(label=labs, effect_measure=measure, lcl=lower, ucl=upper)
    p.labels(effectmeasure='RR')
    p.colors(pointshape="D")
    ax = p.plot(figsize=(7, 3), t_adjuster=0.09, max_value=4, min_value=0.35)
    plt.title("Random Effect Model(Risk Ratio)", loc="right", x=1, y=1.045)
    plt.suptitle("Missing Data Imputation Method", x=-0.1, y=0.98)
    ax.set_xlabel("Favours Control      Favours Haloperidol       ", fontsize=10)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(True)
    ax.spines['left'].set_visible(False)
    plt.savefig("Missing Data Imputation Method", bbox_inches='tight')


if __name__ == '__main__':
    go()