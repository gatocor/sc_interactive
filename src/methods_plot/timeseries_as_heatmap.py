
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

def timeseries_as_heatmap_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.timeseries_as_heatmap(
        config.adata,
        var_names=type_formater(kwargs["var_names"],typing.Collection[str]),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

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

    function = timeseries_as_heatmap_plot
)
