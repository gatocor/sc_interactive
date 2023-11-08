
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

def dpt_groups_pseudotime_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.dpt_groups_pseudotime(
        config.adata,
        color_map=type_formater(kwargs["color_map"],typing.Union[str, matplotlib.colors.Colormap, str]),
        palette=type_formater(kwargs["palette"],typing.Union[typing.Sequence[str], cycler.Cycler, str]),
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

    function = dpt_groups_pseudotime_plot
)