
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

def tracksplot_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.tracksplot(
        config.adata,
        var_names=type_formater(kwargs["var_names"],typing.Union[str, typing.Sequence[str], typing.Mapping[str, typing.Union[str, typing.Sequence[str]]]]),
        groupby=type_formater(kwargs["groupby"],typing.Union[str, typing.Sequence[str]]),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        log=type_formater(kwargs["log"],bool),
        dendrogram=type_formater(kwargs["dendrogram"],typing.Union[bool, str]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        var_group_positions=type_formater(kwargs["var_group_positions"],typing.Optional[typing.Sequence[typing.Tuple[int, int]]]),
        var_group_labels=type_formater(kwargs["var_group_labels"],typing.Optional[typing.Sequence[str]]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        figsize=type_formater(kwargs["figsize"],typing.Optional[typing.Tuple[float, float]]),
    )


    return fig

config.methods_plot["tracksplot"] = dict(
    
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
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['layer'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='figsize', 
        description="typing.Optional[typing.Tuple[float, float]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['figsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = tracksplot_plot,

    docs = sc.pl.umap.__doc__
)
