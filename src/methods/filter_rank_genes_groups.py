
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

filter_rank_genes_groups_args = dict(
    execution = [ARGINPUT,],
    postexecution = [],
    plot = []
)

def filter_rank_genes_groups_f(adata,kwargs):

    sc.tl.filter_rank_genes_groups(
        adata,
    )
        
    return

def filter_rank_genes_groups_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["filter_rank_genes_groups"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_rank_genes_groups_args,

    function = filter_rank_genes_groups_f,

    plot = filter_rank_genes_groups_plot,

)
