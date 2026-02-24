from pathlib import Path
from typing import Literal, Mapping, Union

import matplotlib
import matplotlib.figure
from plotly.graph_objs import Figure


def normalize_filename(filename: str) -> str:
    return filename.replace("/", "_").replace(" ", "_")


PlotlyWriteMode = Literal["cdn"]


def write_plots_to_dir(
    path: Path,
    plots: Mapping[str, Union[Figure, matplotlib.figure.Figure]],
    plotly_js_mode: bool | PlotlyWriteMode = "cdn",
):
    path.mkdir(exist_ok=True, parents=True)
    for name, plot in plots.items():
        if isinstance(plot, Figure):
            target = str(path / normalize_filename(name)) + ".html"
            plot.write_html(target, include_plotlyjs=plotly_js_mode)
            print(f"wrote to {target}")
        elif isinstance(plot, matplotlib.figure.Figure):
            target = str(path / normalize_filename(name)) + ".png"
            plot.savefig(target)
            print(f"wrote to {target}")
        else:
            raise TypeError(f"plot type {type(plot)} not supported")
