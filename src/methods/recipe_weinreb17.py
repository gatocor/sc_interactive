
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

recipe_weinreb17_args = dict(
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
        name='mean_threshold', 
        description="<class 'float'>", 
        visible=dict(function="str(0.01)!=get_node(config.selected)['data']['parameters']['mean_threshold'] or config.show_parameters"),
        properties=dict(value="0.01",type="text")
    ),
    dict(
        input='Input', 
        name='cv_threshold', 
        description="<class 'int'>", 
        visible=dict(function="str(2)!=get_node(config.selected)['data']['parameters']['cv_threshold'] or config.show_parameters"),
        properties=dict(value="2",type="text")
    ),
    dict(
        input='Input', 
        name='n_pcs', 
        description="<class 'int'>", 
        visible=dict(function="str(50)!=get_node(config.selected)['data']['parameters']['n_pcs'] or config.show_parameters"),
        properties=dict(value="50",type="text")
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

def recipe_weinreb17_f(adata,kwargs):

    sc.pp.recipe_weinreb17(
        adata,
        log=type_formater(kwargs["log"],bool),
        mean_threshold=type_formater(kwargs["mean_threshold"],float),
        cv_threshold=type_formater(kwargs["cv_threshold"],int),
        n_pcs=type_formater(kwargs["n_pcs"],int),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def recipe_weinreb17_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["recipe_weinreb17"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = recipe_weinreb17_args,

    function = recipe_weinreb17_f,

    plot = recipe_weinreb17_plot,

)
