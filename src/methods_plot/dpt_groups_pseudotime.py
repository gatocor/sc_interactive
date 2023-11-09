
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

def dpt_groups_pseudotime_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.dpt_groups_pseudotime(
        config.adata,
        color_map=type_formater(kwargs["color_map"],typing.Union[str, matplotlib.colors.Colormap, str]),
        palette=type_formater(kwargs["palette"],typing.Union[typing.Sequence[str], cycler.Cycler, str]),
    )


    return fig

config.methods_plot["dpt_groups_pseudotime"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[str, matplotlib.colors.Colormap, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color_map'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='palette', 
        description="typing.Union[typing.Sequence[str], cycler.Cycler, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['palette'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = dpt_groups_pseudotime_plot
)
