
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

def correlation_matrix_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.correlation_matrix(
        config.adata,
        groupby=type_formater(kwargs["groupby"],str),
        show_correlation_numbers=type_formater(kwargs["show_correlation_numbers"],bool),
        dendrogram=type_formater(kwargs["dendrogram"],typing.Union[bool, str, str]),
        figsize=type_formater(kwargs["figsize"],typing.Optional[typing.Tuple[float, float]]),
        ax=ax,
        vmin=type_formater(kwargs["vmin"],typing.Optional[float]),
        vmax=ax,
        vcenter=type_formater(kwargs["vcenter"],typing.Optional[float]),
        norm=type_formater(kwargs["norm"],typing.Optional[matplotlib.colors.Normalize]),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["correlation_matrix"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='groupby', 
        description="<class 'str'>", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
    ),
    dict(
        input='Input', 
        name='show_correlation_numbers', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['show_correlation_numbers'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='dendrogram', 
        description="typing.Union[bool, str, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['dendrogram'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='figsize', 
        description="typing.Optional[typing.Tuple[float, float]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['figsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vmin', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vmin'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vmax', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vmax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vcenter', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['vcenter'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='norm', 
        description="typing.Optional[matplotlib.colors.Normalize]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['norm'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = correlation_matrix_plot
)
