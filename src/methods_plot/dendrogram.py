
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

def dendrogram_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.dendrogram(
        config.adata,
        groupby=type_formater(kwargs["groupby"],str),
        dendrogram_key=type_formater(kwargs["dendrogram_key"],typing.Optional[str]),
        orientation=type_formater(kwargs["orientation"],typing.Literal['top', 'bottom', 'left', 'right']),
        remove_labels=type_formater(kwargs["remove_labels"],bool),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[str, bool, str]),
        ax=ax,
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png", transparent=True)
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["dendrogram"] = dict(
    
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
        name='dendrogram_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['dendrogram_key'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='orientation', 
        description="typing.Literal['top', 'bottom', 'left', 'right']", 
        visible=dict(function="'top'!=eval(config.active_plot_parameters['orientation']) or config.show_plot"),
        properties=dict(value="'top'",type="text")
    ),
    dict(
        input='Input', 
        name='remove_labels', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['remove_labels'] or config.show_plot"),
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
        description="typing.Union[str, bool, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['save'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = dendrogram_plot
)
