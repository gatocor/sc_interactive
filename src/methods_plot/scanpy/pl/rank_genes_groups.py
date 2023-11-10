
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

def rank_genes_groups_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.rank_genes_groups(
        config.adata,
        groups=type_formater(kwargs["groups"],typing.Union[str, typing.Sequence[str]]),
        n_genes=type_formater(kwargs["n_genes"],int),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        key=type_formater(kwargs["key"],typing.Optional[str]),
        fontsize=type_formater(kwargs["fontsize"],int),
        ncols=type_formater(kwargs["ncols"],int),
        sharey=type_formater(kwargs["sharey"],bool),
        ax=ax,
    )


    return fig

config.methods_plot["rank_genes_groups"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='groups', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groups'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_genes', 
        description="<class 'int'>", 
        visible=dict(function="str(20)!=config.active_plot_parameters['n_genes'] or config.show_plot"),
        properties=dict(value="20",type="text")
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
        name='key', 
        description="typing.Optional[str]", 
        visible=dict(function="'rank_genes_groups'!=eval(config.active_plot_parameters['key']) or config.show_plot"),
        properties=dict(value="'rank_genes_groups'",type="text")
    ),
    dict(
        input='Input', 
        name='fontsize', 
        description="<class 'int'>", 
        visible=dict(function="str(8)!=config.active_plot_parameters['fontsize'] or config.show_plot"),
        properties=dict(value="8",type="text")
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
        name='sharey', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['sharey'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = rank_genes_groups_plot,

    docs = sc.pl.umap.__doc__
)
