
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

def highest_expr_genes_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.highest_expr_genes(
        config.adata,
        n_top=type_formater(kwargs["n_top"],int),
        ax=ax,
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        log=type_formater(kwargs["log"],bool),
    )


    return fig

config.methods_plot["highest_expr_genes"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='n_top', 
        description="<class 'int'>", 
        visible=dict(function="str(30)!=config.active_plot_parameters['n_top'] or config.show_plot"),
        properties=dict(value="30",type="text")
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
        name='gene_symbols', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['gene_symbols'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),],

    function = highest_expr_genes_plot
)
