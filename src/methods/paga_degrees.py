
import numpy
from numpy import inf
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scanpy as sc
import louvain
import scipy
import leidenalg
import plotly.tools as tls
import cycler
import matplotlib      # pip install matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO
from general import *

paga_degrees_args = dict(
    execution = [ARGINPUT,],
    postexecution = [],
    plot = []
)

def paga_degrees_f(adata,kwargs):

    sc.tl.paga_degrees(
        adata,
    )
        
    return

def paga_degrees_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["paga_degrees"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = paga_degrees_args,

    function = paga_degrees_f,

    plot = paga_degrees_plot,

)
