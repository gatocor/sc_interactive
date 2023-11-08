
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

combat_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='key', 
        description="<class 'str'>", 
        visible=dict(function="'batch'!=eval(config.active_node_parameters['key']) or config.show_parameters"),
        properties=dict(value="'batch'",type="text")
    ),
    dict(
        input='Input', 
        name='covariates', 
        description="typing.Optional[typing.Collection[str]]", 
        visible=dict(function="str(None)!=config.active_node_parameters['covariates'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='inplace', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['inplace'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),],
    postexecution = [],
    plot = []
)

def combat_f(adata,kwargs):

    sc.pp.combat(
        adata,
        key=type_formater(kwargs["key"],str),
        covariates=type_formater(kwargs["covariates"],typing.Optional[typing.Collection[str]]),
        inplace=type_formater(kwargs["inplace"],bool),
    )
        
    return

def combat_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["combat"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = combat_args,

    function = combat_f,

    plot = combat_plot,

)
