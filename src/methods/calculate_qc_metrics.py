
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

calculate_qc_metrics_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='expr_type', 
        description="<class 'str'>", 
        visible=dict(function="'counts'!=eval(config.active_node_parameters['expr_type']) or config.show_parameters"),
        properties=dict(value="'counts'",type="text")
    ),
    dict(
        input='Input', 
        name='var_type', 
        description="<class 'str'>", 
        visible=dict(function="'genes'!=eval(config.active_node_parameters['var_type']) or config.show_parameters"),
        properties=dict(value="'genes'",type="text")
    ),
    dict(
        input='Input', 
        name='qc_vars', 
        description="typing.Collection[str]", 
        visible=dict(function="str(())!=config.active_node_parameters['qc_vars'] or config.show_parameters"),
        properties=dict(value="()",type="text")
    ),
    dict(
        input='Input', 
        name='percent_top', 
        description="typing.Optional[typing.Collection[int]]", 
        visible=dict(function="str((50, 100, 200, 500))!=config.active_node_parameters['percent_top'] or config.show_parameters"),
        properties=dict(value="(50, 100, 200, 500)",type="text")
    ),
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['layer'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_raw', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['use_raw'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='inplace', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['inplace'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='log1p', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['log1p'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='parallel', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_node_parameters['parallel'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),],
    postexecution = [],
    plot = []
)

def calculate_qc_metrics_f(adata,kwargs):

    sc.pp.calculate_qc_metrics(
        adata,
        expr_type=type_formater(kwargs["expr_type"],str),
        var_type=type_formater(kwargs["var_type"],str),
        qc_vars=type_formater(kwargs["qc_vars"],typing.Collection[str]),
        percent_top=type_formater(kwargs["percent_top"],typing.Optional[typing.Collection[int]]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        use_raw=type_formater(kwargs["use_raw"],bool),
        inplace=type_formater(kwargs["inplace"],bool),
        log1p=type_formater(kwargs["log1p"],bool),
        parallel=type_formater(kwargs["parallel"],typing.Optional[bool]),
    )
        
    return

def calculate_qc_metrics_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["calculate_qc_metrics"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = calculate_qc_metrics_args,

    function = calculate_qc_metrics_f,

    plot = calculate_qc_metrics_plot,

)
