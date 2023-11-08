
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

sqrt_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='chunked', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['chunked'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='chunk_size', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['chunk_size'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),],
    postexecution = [],
    plot = []
)

def sqrt_f(adata,kwargs):

    sc.pp.sqrt(
        adata,
        copy=type_formater(kwargs["copy"],bool),
        chunked=type_formater(kwargs["chunked"],bool),
        chunk_size=type_formater(kwargs["chunk_size"],typing.Optional[int]),
    )
        
    return

def sqrt_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["sqrt"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = sqrt_args,

    function = sqrt_f,

    plot = sqrt_plot,

)
