
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

def highly_variable_genes_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.highly_variable_genes(
        config.adata,
        log=type_formater(kwargs["log"],bool),
        highly_variable_genes=type_formater(kwargs["highly_variable_genes"],bool),
    )


    return fig

config.methods_plot["highly_variable_genes"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='highly_variable_genes', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['highly_variable_genes'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),],

    function = highly_variable_genes_plot
)
