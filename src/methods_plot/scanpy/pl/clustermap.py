
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

def clustermap_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.clustermap(
        config.adata,
        obs_keys=type_formater(kwargs["obs_keys"],str),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
    )


    return fig

config.methods_plot["clustermap"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='obs_keys', 
        description="<class 'str'>", 
        visible=dict(function="str(None)!=config.active_plot_parameters['obs_keys'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_raw', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['use_raw'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = clustermap_plot
)
