
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

def pca_variance_ratio_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.pca_variance_ratio(
        config.adata,
        n_pcs=type_formater(kwargs["n_pcs"],int),
        log=type_formater(kwargs["log"],bool),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png", transparent=True)
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["pca_variance_ratio"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='n_pcs', 
        description="<class 'int'>", 
        visible=dict(function="str(30)!=config.active_plot_parameters['n_pcs'] or config.show_plot"),
        properties=dict(value="30",type="text")
    ),
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),],

    function = pca_variance_ratio_plot
)
