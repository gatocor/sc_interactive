
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

subsample_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='fraction', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['fraction'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_obs', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['n_obs'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='random_state', 
        description="typing.Union[str, int, numpy.random.mtrand.RandomState]", 
        visible=dict(function="str(0)!=get_node(config.selected)['data']['parameters']['random_state'] or config.show_parameters"),
        properties=dict(value="0",type="text")
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

def subsample_f(adata,kwargs):

    sc.pp.subsample(
        adata,
        fraction=type_formater(kwargs["fraction"],typing.Optional[float]),
        n_obs=type_formater(kwargs["n_obs"],typing.Optional[int]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def subsample_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["subsample"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = subsample_args,

    function = subsample_f,

    plot = subsample_plot,

)
