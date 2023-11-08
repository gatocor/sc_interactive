
import numpy
from numpy import inf
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scipy
import plotly.tools as tls
import cycler
import matplotlib      # pip install matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO

from general import *

def pca_scatter_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.pca_scatter(
        config.adata,
        color=type_formater(kwargs["color"],typing.Union[str, typing.Sequence[str], str]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        sort_order=type_formater(kwargs["sort_order"],bool),
        edges=type_formater(kwargs["edges"],bool),
        edges_width=type_formater(kwargs["edges_width"],float),
        edges_color=type_formater(kwargs["edges_color"],typing.Union[str, typing.Sequence[float], typing.Sequence[str]]),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        arrows=type_formater(kwargs["arrows"],bool),
        arrows_kwds=type_formater(kwargs["arrows_kwds"],typing.Optional[typing.Mapping[str, typing.Any]]),
        groups=type_formater(kwargs["groups"],typing.Optional[str]),
        components=type_formater(kwargs["components"],typing.Union[str, typing.Sequence[str]]),
        dimensions=type_formater(kwargs["dimensions"],typing.Union[typing.Tuple[int, int], typing.Sequence[typing.Tuple[int, int]], str]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        projection=type_formater(kwargs["projection"],typing.Literal['2d', '3d']),
        scale_factor=type_formater(kwargs["scale_factor"],typing.Optional[float]),
        color_map=type_formater(kwargs["color_map"],typing.Union[matplotlib.colors.Colormap, str, str]),
        cmap=type_formater(kwargs["cmap"],typing.Union[matplotlib.colors.Colormap, str, str]),
        palette=type_formater(kwargs["palette"],typing.Union[str, typing.Sequence[str], cycler.Cycler, str]),
        na_color=type_formater(kwargs["na_color"],typing.Union[str, typing.Tuple[float, ...]]),
        na_in_legend=type_formater(kwargs["na_in_legend"],bool),
        size=type_formater(kwargs["size"],typing.Union[float, typing.Sequence[float], str]),
        frameon=type_formater(kwargs["frameon"],typing.Optional[bool]),
        legend_fontsize=type_formater(kwargs["legend_fontsize"],typing.Union[int, float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], str]),
        legend_fontweight=type_formater(kwargs["legend_fontweight"],typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]),
        legend_loc=type_formater(kwargs["legend_loc"],str),
        legend_fontoutline=type_formater(kwargs["legend_fontoutline"],typing.Optional[int]),
        colorbar_loc=type_formater(kwargs["colorbar_loc"],typing.Optional[str]),
        vmax=type_formater(kwargs["vmax"],typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]),
        vmin=type_formater(kwargs["vmin"],typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]),
        vcenter=type_formater(kwargs["vcenter"],typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]),
        norm=type_formater(kwargs["norm"],typing.Union[matplotlib.colors.Normalize, typing.Sequence[matplotlib.colors.Normalize], str]),
        add_outline=type_formater(kwargs["add_outline"],typing.Optional[bool]),
        outline_width=type_formater(kwargs["outline_width"],typing.Tuple[float, float]),
        outline_color=type_formater(kwargs["outline_color"],typing.Tuple[str, str]),
        ncols=type_formater(kwargs["ncols"],int),
        hspace=type_formater(kwargs["hspace"],float),
        wspace=type_formater(kwargs["wspace"],typing.Optional[float]),
        title=type_formater(kwargs["title"],typing.Union[str, typing.Sequence[str], str]),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
        ax=type_formater(kwargs["ax"],typing.Optional[matplotlib.axes._axes.Axes]),
        return_fig=True,
        annotate_var_explained=type_formater(kwargs["annotate_var_explained"],bool),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["pca_scatter"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='min_dist', 
        description="<class 'float'>", 
        visible=dict(function="str(0.5)!=config.active_node_parameters['min_dist'] or config.show_parameters"),
        properties=dict(value="0.5",type="text")
    ),
    dict(
        input='Input', 
        name='spread', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['spread'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='n_components', 
        description="<class 'int'>", 
        visible=dict(function="str(2)!=config.active_node_parameters['n_components'] or config.show_parameters"),
        properties=dict(value="2",type="text")
    ),
    dict(
        input='Input', 
        name='maxiter', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['maxiter'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='alpha', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['alpha'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='gamma', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['gamma'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='negative_sample_rate', 
        description="<class 'int'>", 
        visible=dict(function="str(5)!=config.active_node_parameters['negative_sample_rate'] or config.show_parameters"),
        properties=dict(value="5",type="text")
    ),
    dict(
        input='Input', 
        name='init_pos', 
        description="typing.Union[typing.Literal['paga', 'spectral', 'random'], numpy.ndarray, str]", 
        visible=dict(function="'spectral'!=eval(config.active_node_parameters['init_pos']) or config.show_parameters"),
        properties=dict(value="'spectral'",type="text")
    ),
    dict(
        input='Input', 
        name='random_state', 
        description="typing.Union[str, int, numpy.random.mtrand.RandomState]", 
        visible=dict(function="str(0)!=config.active_node_parameters['random_state'] or config.show_parameters"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='a', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['a'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='b', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['b'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='method', 
        description="typing.Literal['umap', 'rapids']", 
        visible=dict(function="'umap'!=eval(config.active_node_parameters['method']) or config.show_parameters"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['neighbors_key'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),],

    function = pca_scatter_plot
)
