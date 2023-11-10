
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

def filter_genes_dispersion_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.filter_genes_dispersion(
        config.adata,
        log=type_formater(kwargs["log"],bool),
    )


    return fig

config.methods_plot["filter_genes_dispersion"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),],

    function = filter_genes_dispersion_plot,

    docs = sc.pl.umap.__doc__
)
