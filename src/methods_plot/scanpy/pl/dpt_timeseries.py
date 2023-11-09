
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

def dpt_timeseries_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.dpt_timeseries(
        config.adata,
        color_map=type_formater(kwargs["color_map"],typing.Union[str, matplotlib.colors.Colormap]),
        as_heatmap=type_formater(kwargs["as_heatmap"],bool),
    )


    return fig

config.methods_plot["dpt_timeseries"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[str, matplotlib.colors.Colormap]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color_map'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='as_heatmap', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['as_heatmap'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),],

    function = dpt_timeseries_plot
)