
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

def sim_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.sim(
        config.adata,
        tmax_realization=type_formater(kwargs["tmax_realization"],typing.Optional[int]),
        as_heatmap=type_formater(kwargs["as_heatmap"],bool),
        shuffle=type_formater(kwargs["shuffle"],bool),
    )


    return fig

config.methods_plot["sim"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='tmax_realization', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['tmax_realization'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='as_heatmap', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['as_heatmap'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='shuffle', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['shuffle'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),],

    function = sim_plot,

    docs = sc.pl.umap.__doc__
)
