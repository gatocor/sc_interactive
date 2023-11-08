
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

filter_genes_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='min_counts', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['min_counts'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='min_cells', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['min_cells'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_counts', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['max_counts'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_cells', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['max_cells'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='inplace', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['inplace'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),],
    postexecution = [],
    plot = [
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='show', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_node_parameters['show'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='save', 
        description="typing.Union[bool, str, str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['save'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),]
)

def filter_genes_f(adata,kwargs):

    sc.pp.filter_genes(
        adata,
        min_counts=type_formater(kwargs["min_counts"],typing.Optional[int]),
        min_cells=type_formater(kwargs["min_cells"],typing.Optional[int]),
        max_counts=type_formater(kwargs["max_counts"],typing.Optional[int]),
        max_cells=type_formater(kwargs["max_cells"],typing.Optional[int]),
        inplace=type_formater(kwargs["inplace"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def filter_genes_plot():

    kwargs = get_node(config.selected)['data']['plot']
    
    fig = sc.pl.filter_genes_dispersion(
        config.adata,
        log=type_formater(kwargs["log"],bool),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
    )


    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["filter_genes"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_genes_args,

    function = filter_genes_f,

    plot = filter_genes_plot,

)
