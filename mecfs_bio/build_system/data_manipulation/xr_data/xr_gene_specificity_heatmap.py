import plotly.express as px
import structlog
import xarray as xr
from attrs import frozen
from plotly.graph_objs import Figure

from mecfs_bio.constants.xr_constants import XR_SPECIFICITY_MATRIX, XR_TISSUE_DIMENSION


@frozen
class GeneNormalizedPlotMode:
    max_multiple: float


XRHeatmapPlotMode = GeneNormalizedPlotMode


@frozen
class XRHeatmapPlotSpec:
    plot_mode: XRHeatmapPlotMode
    color_scale: str | list[str]
    y_array: str | None = None


logger = structlog.get_logger()


def xr_make_gene_spec_heatmap(ds: xr.Dataset, plot_spec: XRHeatmapPlotSpec) -> Figure:
    """
    Generate a heatmap illustrating a gene/tissue specificity matrix45G
    """
    ds = ds.copy()
    plot_mode = plot_spec.plot_mode
    logger.debug(f"plotting gene/tissue dataset of dimension {ds.sizes}")
    plot_options: dict = {}
    plot_options["aspect"] = "auto"
    if isinstance(plot_mode, GeneNormalizedPlotMode):
        ds[XR_SPECIFICITY_MATRIX] = ds[XR_SPECIFICITY_MATRIX] / ds[
            XR_SPECIFICITY_MATRIX
        ].sum(dim=XR_TISSUE_DIMENSION)
        num_tissues = ds.sizes[XR_TISSUE_DIMENSION]
        expected_frac = 1 / num_tissues
        zmin = expected_frac - plot_mode.max_multiple * expected_frac
        zmax = expected_frac + plot_mode.max_multiple * expected_frac
        plot_options["zmin"] = zmin
        plot_options["zmax"] = zmax
        plot_options["color_continuous_scale"] = plot_spec.color_scale  # "Tropic"
    if plot_spec.y_array is not None:
        plot_options["y"] = ds[plot_spec.y_array]
    fig = px.imshow(ds[XR_SPECIFICITY_MATRIX], **plot_options)
    return fig
