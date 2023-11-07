
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

filter_cells_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='min_counts', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['min_counts']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='min_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['min_genes']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_counts', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['max_counts']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['max_genes']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='inplace', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['inplace']"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['copy']"),
        properties=dict(value="False",type="text")
    ),],
    postexecution = [],
    plot = []
)

def filter_cells_f(adata,kwargs):

    sc.pp.filter_cells(
        adata,
        min_counts=type_formater(kwargs["min_counts"],typing.Optional[int]),
        min_genes=type_formater(kwargs["min_genes"],typing.Optional[int]),
        max_counts=type_formater(kwargs["max_counts"],typing.Optional[int]),
        max_genes=type_formater(kwargs["max_genes"],typing.Optional[int]),
        inplace=type_formater(kwargs["inplace"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def filter_cells_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["filter_cells"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_cells_args,

    function = filter_cells_f,

    plot = filter_cells_plot,

)
