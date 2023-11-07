
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

recipe_zheng17_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='n_top_genes', 
        description="<class 'int'>", 
        visible=dict(function="str(1000)!=get_node(config.selected)['data']['parameters']['n_top_genes'] or config.show_parameters"),
        properties=dict(value="1000",type="text")
    ),
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

def recipe_zheng17_f(adata,kwargs):

    sc.pp.recipe_zheng17(
        adata,
        n_top_genes=type_formater(kwargs["n_top_genes"],int),
        log=type_formater(kwargs["log"],bool),
        plot=type_formater(kwargs["plot"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def recipe_zheng17_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["recipe_zheng17"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = recipe_zheng17_args,

    function = recipe_zheng17_f,

    plot = recipe_zheng17_plot,

)
