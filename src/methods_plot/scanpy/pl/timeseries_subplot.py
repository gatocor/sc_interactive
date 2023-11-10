
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

def timeseries_subplot_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.timeseries_subplot(
        config.adata,
        palette=type_formater(kwargs["palette"],typing.Union[typing.Sequence[str], cycler.Cycler, str]),
        ax=ax,
    )


    return fig

config.methods_plot["timeseries_subplot"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='palette', 
        description="typing.Union[typing.Sequence[str], cycler.Cycler, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['palette'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = timeseries_subplot_plot,

    docs = sc.pl.umap.__doc__
)
