
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

def paga_adjacency_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.paga_adjacency(
        config.adata,
    )


    return fig

config.methods_plot["paga_adjacency"] = dict(
    
    args = [],

    function = paga_adjacency_plot
)
