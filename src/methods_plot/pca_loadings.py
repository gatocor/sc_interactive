
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

def pca_loadings_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.pca_loadings(
        config.adata,
        components=type_formater(kwargs["components"],typing.Union[str, typing.Sequence[int], str]),
        include_lowest=type_formater(kwargs["include_lowest"],bool),
        n_points=type_formater(kwargs["n_points"],typing.Optional[int]),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[str, bool, str]),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["pca_loadings"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='components', 
        description="typing.Union[str, typing.Sequence[int], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['components'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='include_lowest', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['include_lowest'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='n_points', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['n_points'] or config.show_plot"),
        properties=dict(value="None",type="text")
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
    ),],

    function = pca_loadings_plot
)
