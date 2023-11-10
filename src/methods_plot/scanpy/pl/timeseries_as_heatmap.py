
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

def timeseries_as_heatmap_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.timeseries_as_heatmap(
        config.adata,
        var_names=type_formater(kwargs["var_names"],typing.Collection[str]),
    )


    return fig

config.methods_plot["timeseries_as_heatmap"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='var_names', 
        description="typing.Collection[str]", 
        visible=dict(function="str(())!=config.active_plot_parameters['var_names'] or config.show_plot"),
        properties=dict(value="()",type="text")
    ),],

    function = timeseries_as_heatmap_plot,

    docs = sc.pl.umap.__doc__
)
