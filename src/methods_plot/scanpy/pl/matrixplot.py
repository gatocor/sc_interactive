
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

def matrixplot_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.matrixplot(
        config.adata,
        var_names=type_formater(kwargs["var_names"],typing.Union[str, typing.Sequence[str], typing.Mapping[str, typing.Union[str, typing.Sequence[str]]]]),
        groupby=type_formater(kwargs["groupby"],typing.Union[str, typing.Sequence[str]]),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        log=type_formater(kwargs["log"],bool),
        num_categories=type_formater(kwargs["num_categories"],int),
        figsize=type_formater(kwargs["figsize"],typing.Optional[typing.Tuple[float, float]]),
        dendrogram=type_formater(kwargs["dendrogram"],typing.Union[bool, str]),
        title=type_formater(kwargs["title"],typing.Optional[str]),
        cmap=type_formater(kwargs["cmap"],typing.Optional[str]),
        colorbar_title=type_formater(kwargs["colorbar_title"],typing.Optional[str]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        var_group_positions=type_formater(kwargs["var_group_positions"],typing.Optional[typing.Sequence[typing.Tuple[int, int]]]),
        var_group_labels=type_formater(kwargs["var_group_labels"],typing.Optional[typing.Sequence[str]]),
        var_group_rotation=type_formater(kwargs["var_group_rotation"],typing.Optional[float]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        standard_scale=type_formater(kwargs["standard_scale"],typing.Literal['var', 'group']),
        values_df=type_formater(kwargs["values_df"],typing.Optional[pandas.core.frame.DataFrame]),
        swap_axes=type_formater(kwargs["swap_axes"],bool),
        ax=ax,
        vmin=type_formater(kwargs["vmin"],typing.Optional[float]),
        vmax=ax,
        vcenter=type_formater(kwargs["vcenter"],typing.Optional[float]),
        norm=type_formater(kwargs["norm"],typing.Optional[matplotlib.colors.Normalize]),
    )


    return fig

config.methods_plot["matrixplot"] = dict(
    
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
        name='use_raw', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['use_raw'] or config.show_plot"),
        properties=dict(value="None",type="text")
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
        name='num_categories', 
        description="<class 'int'>", 
        visible=dict(function="str(7)!=config.active_plot_parameters['num_categories'] or config.show_plot"),
        properties=dict(value="7",type="text")
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
        name='title', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['title'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cmap', 
        description="typing.Optional[str]", 
        visible=dict(function="'viridis'!=eval(config.active_plot_parameters['cmap']) or config.show_plot"),
        properties=dict(value="'viridis'",type="text")
    ),
    dict(
        input='Input', 
        name='colorbar_title', 
        description="typing.Optional[str]", 
        visible=dict(function="'Mean expression in group'!=eval(config.active_plot_parameters['colorbar_title']) or config.show_plot"),
        properties=dict(value="'Mean expression in group'",type="text")
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
        name='standard_scale', 
        description="typing.Literal['var', 'group']", 
        visible=dict(function="str(None)!=config.active_plot_parameters['standard_scale'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='values_df', 
        description="typing.Optional[pandas.core.frame.DataFrame]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['values_df'] or config.show_plot"),
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

    function = matrixplot_plot,

    docs = sc.pl.umap.__doc__
)
