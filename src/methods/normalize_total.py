
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

normalize_total_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='target_sum', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['target_sum'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='exclude_highly_expressed', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['exclude_highly_expressed'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='max_fraction', 
        description="<class 'float'>", 
        visible=dict(function="str(0.05)!=get_node(config.selected)['data']['parameters']['max_fraction'] or config.show_parameters"),
        properties=dict(value="0.05",type="text")
    ),
    dict(
        input='Input', 
        name='key_added', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['key_added'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['layer'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layers', 
        description="typing.Union[typing.Literal['all'], typing.Iterable[str]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['layers'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layer_norm', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['layer_norm'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='inplace', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['inplace'] or config.show_parameters"),
        properties=dict(value="True",type="text")
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

def normalize_total_f(adata,kwargs):

    sc.pp.normalize_total(
        adata,
        target_sum=type_formater(kwargs["target_sum"],typing.Optional[float]),
        exclude_highly_expressed=type_formater(kwargs["exclude_highly_expressed"],bool),
        max_fraction=type_formater(kwargs["max_fraction"],float),
        key_added=type_formater(kwargs["key_added"],typing.Optional[str]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        layers=type_formater(kwargs["layers"],typing.Union[typing.Literal['all'], typing.Iterable[str]]),
        layer_norm=type_formater(kwargs["layer_norm"],typing.Optional[str]),
        inplace=type_formater(kwargs["inplace"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def normalize_total_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["normalize_total"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = normalize_total_args,

    function = normalize_total_f,

    plot = normalize_total_plot,

)
