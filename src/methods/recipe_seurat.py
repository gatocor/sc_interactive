
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

recipe_seurat_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['log'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='plot', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['plot'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),],
    postexecution = [],
    plot = []
)

def recipe_seurat_f(adata,kwargs):

    sc.pp.recipe_seurat(
        adata,
        log=type_formater(kwargs["log"],bool),
        plot=type_formater(kwargs["plot"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def recipe_seurat_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["recipe_seurat"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = recipe_seurat_args,

    function = recipe_seurat_f,

    plot = recipe_seurat_plot,

)