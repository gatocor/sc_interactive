
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
import base64
from io import BytesIO

from general import *

def sim_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.sim(
        config.adata,
        tmax_realization=type_formater(kwargs["tmax_realization"],typing.Optional[int]),
        as_heatmap=type_formater(kwargs["as_heatmap"],bool),
        shuffle=type_formater(kwargs["shuffle"],bool),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

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
    ),
    dict(
        input='Input', 
        name='show', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['show'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='save', 
        description="typing.Union[bool, str, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['save'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = sim_plot
)