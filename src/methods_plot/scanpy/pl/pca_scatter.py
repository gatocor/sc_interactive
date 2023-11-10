
import numpy
from numpy import inf
import scanpy
import pandas
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scipy
import plotly.tools as tls
import cycler
import matplotlib      # pip install matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from general import *

def pca_scatter_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.pca_scatter(
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
        vmax=ax,
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
        ax=ax,
        return_fig=False,
        annotate_var_explained=type_formater(kwargs["annotate_var_explained"],bool),
    )


    return fig

config.methods_plot["pca_scatter"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='color', 
        description="typing.Union[str, typing.Sequence[str], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='gene_symbols', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['gene_symbols'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_raw', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['use_raw'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='sort_order', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['sort_order'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='edges', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['edges'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='edges_width', 
        description="<class 'float'>", 
        visible=dict(function="str(0.1)!=config.active_plot_parameters['edges_width'] or config.show_plot"),
        properties=dict(value="0.1",type="text")
    ),
    dict(
        input='Input', 
        name='edges_color', 
        description="typing.Union[str, typing.Sequence[float], typing.Sequence[str]]", 
        visible=dict(function="'grey'!=eval(config.active_plot_parameters['edges_color']) or config.show_plot"),
        properties=dict(value="'grey'",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['neighbors_key'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='arrows', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['arrows'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='arrows_kwds', 
        description="typing.Optional[typing.Mapping[str, typing.Any]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['arrows_kwds'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='groups', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groups'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='components', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['components'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='dimensions', 
        description="typing.Union[typing.Tuple[int, int], typing.Sequence[typing.Tuple[int, int]], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['dimensions'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['layer'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='projection', 
        description="typing.Literal['2d', '3d']", 
        visible=dict(function="'2d'!=eval(config.active_plot_parameters['projection']) or config.show_plot"),
        properties=dict(value="'2d'",type="text")
    ),
    dict(
        input='Input', 
        name='scale_factor', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['scale_factor'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[matplotlib.colors.Colormap, str, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color_map'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cmap', 
        description="typing.Union[matplotlib.colors.Colormap, str, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['cmap'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='palette', 
        description="typing.Union[str, typing.Sequence[str], cycler.Cycler, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['palette'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='na_color', 
        description="typing.Union[str, typing.Tuple[float, ...]]", 
        visible=dict(function="'lightgray'!=eval(config.active_plot_parameters['na_color']) or config.show_plot"),
        properties=dict(value="'lightgray'",type="text")
    ),
    dict(
        input='Input', 
        name='na_in_legend', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['na_in_legend'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='size', 
        description="typing.Union[float, typing.Sequence[float], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['size'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='frameon', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['frameon'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontsize', 
        description="typing.Union[int, float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['legend_fontsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontweight', 
        description="typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]", 
        visible=dict(function="'bold'!=eval(config.active_plot_parameters['legend_fontweight']) or config.show_plot"),
        properties=dict(value="'bold'",type="text")
    ),
    dict(
        input='Input', 
        name='legend_loc', 
        description="<class 'str'>", 
        visible=dict(function="'right margin'!=eval(config.active_plot_parameters['legend_loc']) or config.show_plot"),
        properties=dict(value="'right margin'",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontoutline', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['legend_fontoutline'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='colorbar_loc', 
        description="typing.Optional[str]", 
        visible=dict(function="'right'!=eval(config.active_plot_parameters['colorbar_loc']) or config.show_plot"),
        properties=dict(value="'right'",type="text")
    ),
    dict(
        input='Input', 
        name='vmax', 
        description="typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vmax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vmin', 
        description="typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vmin'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vcenter', 
        description="typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vcenter'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='norm', 
        description="typing.Union[matplotlib.colors.Normalize, typing.Sequence[matplotlib.colors.Normalize], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['norm'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='add_outline', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(False)!=config.active_plot_parameters['add_outline'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='outline_width', 
        description="typing.Tuple[float, float]", 
        visible=dict(function="str((0.3, 0.05))!=config.active_plot_parameters['outline_width'] or config.show_plot"),
        properties=dict(value="(0.3, 0.05)",type="text")
    ),
    dict(
        input='Input', 
        name='outline_color', 
        description="typing.Tuple[str, str]", 
        visible=dict(function="str(('black', 'white'))!=config.active_plot_parameters['outline_color'] or config.show_plot"),
        properties=dict(value="('black', 'white')",type="text")
    ),
    dict(
        input='Input', 
        name='ncols', 
        description="<class 'int'>", 
        visible=dict(function="str(4)!=config.active_plot_parameters['ncols'] or config.show_plot"),
        properties=dict(value="4",type="text")
    ),
    dict(
        input='Input', 
        name='hspace', 
        description="<class 'float'>", 
        visible=dict(function="str(0.25)!=config.active_plot_parameters['hspace'] or config.show_plot"),
        properties=dict(value="0.25",type="text")
    ),
    dict(
        input='Input', 
        name='wspace', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['wspace'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='title', 
        description="typing.Union[str, typing.Sequence[str], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['title'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='show', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['show'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='save', 
        description="typing.Union[bool, str, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['save'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='return_fig', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['return_fig'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='annotate_var_explained', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['annotate_var_explained'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),],

    function = pca_scatter_plot,

    docs = sc.pl.umap.__doc__
)
