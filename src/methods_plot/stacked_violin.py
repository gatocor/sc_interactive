
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

def stacked_violin_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.stacked_violin(
        config.adata,
        var_names=type_formater(kwargs["var_names"],typing.Union[str, typing.Sequence[str], typing.Mapping[str, typing.Union[str, typing.Sequence[str]]]]),
        groupby=type_formater(kwargs["groupby"],typing.Union[str, typing.Sequence[str]]),
        log=type_formater(kwargs["log"],bool),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        num_categories=type_formater(kwargs["num_categories"],int),
        title=type_formater(kwargs["title"],typing.Optional[str]),
        colorbar_title=type_formater(kwargs["colorbar_title"],typing.Optional[str]),
        figsize=type_formater(kwargs["figsize"],typing.Optional[typing.Tuple[float, float]]),
        dendrogram=type_formater(kwargs["dendrogram"],typing.Union[bool, str]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        var_group_positions=type_formater(kwargs["var_group_positions"],typing.Optional[typing.Sequence[typing.Tuple[int, int]]]),
        var_group_labels=type_formater(kwargs["var_group_labels"],typing.Optional[typing.Sequence[str]]),
        standard_scale=type_formater(kwargs["standard_scale"],typing.Optional[typing.Literal['var', 'obs']]),
        var_group_rotation=type_formater(kwargs["var_group_rotation"],typing.Optional[float]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        stripplot=type_formater(kwargs["stripplot"],bool),
        jitter=type_formater(kwargs["jitter"],typing.Union[float, bool]),
        size=type_formater(kwargs["size"],int),
        scale=type_formater(kwargs["scale"],typing.Literal['area', 'count', 'width']),
        yticklabels=type_formater(kwargs["yticklabels"],typing.Optional[bool]),
        order=type_formater(kwargs["order"],typing.Optional[typing.Sequence[str]]),
        swap_axes=type_formater(kwargs["swap_axes"],bool),
        row_palette=type_formater(kwargs["row_palette"],typing.Optional[str]),
        cmap=type_formater(kwargs["cmap"],typing.Optional[str]),
        ax=ax,
        vmin=type_formater(kwargs["vmin"],typing.Optional[float]),
        vmax=ax,
        vcenter=type_formater(kwargs["vcenter"],typing.Optional[float]),
        norm=type_formater(kwargs["norm"],typing.Optional[matplotlib.colors.Normalize]),
    )


    return fig

config.methods_plot["stacked_violin"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='var_names', 
        description="typing.Union[str, typing.Sequence[str], typing.Mapping[str, typing.Union[str, typing.Sequence[str]]]]", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
    ),
    dict(
        input='Input', 
        name='groupby', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
    ),
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
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
        name='num_categories', 
        description="<class 'int'>", 
        visible=dict(function="str(7)!=config.active_plot_parameters['num_categories'] or config.show_plot"),
        properties=dict(value="7",type="text")
    ),
    dict(
        input='Input', 
        name='title', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['title'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='colorbar_title', 
        description="typing.Optional[str]", 
        visible=dict(function="'Median expression in group'!=eval(config.active_plot_parameters['colorbar_title']) or config.show_plot"),
        properties=dict(value="'Median expression in group'",type="text")
    ),
    dict(
        input='Input', 
        name='figsize', 
        description="typing.Optional[typing.Tuple[float, float]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['figsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='dendrogram', 
        description="typing.Union[bool, str]", 
        visible=dict(function="str(False)!=config.active_plot_parameters['dendrogram'] or config.show_plot"),
        properties=dict(value="False",type="text")
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
        name='var_group_positions', 
        description="typing.Optional[typing.Sequence[typing.Tuple[int, int]]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['var_group_positions'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='var_group_labels', 
        description="typing.Optional[typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['var_group_labels'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='standard_scale', 
        description="typing.Optional[typing.Literal['var', 'obs']]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['standard_scale'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='var_group_rotation', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['var_group_rotation'] or config.show_plot"),
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
        name='stripplot', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['stripplot'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='jitter', 
        description="typing.Union[float, bool]", 
        visible=dict(function="str(False)!=config.active_plot_parameters['jitter'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='size', 
        description="<class 'int'>", 
        visible=dict(function="str(1)!=config.active_plot_parameters['size'] or config.show_plot"),
        properties=dict(value="1",type="text")
    ),
    dict(
        input='Input', 
        name='scale', 
        description="typing.Literal['area', 'count', 'width']", 
        visible=dict(function="'width'!=eval(config.active_plot_parameters['scale']) or config.show_plot"),
        properties=dict(value="'width'",type="text")
    ),
    dict(
        input='Input', 
        name='yticklabels', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(False)!=config.active_plot_parameters['yticklabels'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='order', 
        description="typing.Optional[typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['order'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='swap_axes', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['swap_axes'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='row_palette', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['row_palette'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cmap', 
        description="typing.Optional[str]", 
        visible=dict(function="'Blues'!=eval(config.active_plot_parameters['cmap']) or config.show_plot"),
        properties=dict(value="'Blues'",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[scanpy.plotting._utils._AxesSubplot]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vmin', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vmin'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vmax', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vmax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vcenter', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vcenter'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='norm', 
        description="typing.Optional[matplotlib.colors.Normalize]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['norm'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = stacked_violin_plot
)
