
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

filter_genes_dispersion_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='flavor', 
        description="typing.Literal['seurat', 'cell_ranger']", 
        visible=dict(function="'seurat'!=eval(get_node(config.selected)['data']['parameters']['flavor']) or config.show_parameters"),
        properties=dict(value="'seurat'",type="text")
    ),
    dict(
        input='Input', 
        name='min_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['min_disp'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['max_disp'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='min_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['min_mean'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['max_mean'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_bins', 
        description="<class 'int'>", 
        visible=dict(function="str(20)!=get_node(config.selected)['data']['parameters']['n_bins'] or config.show_parameters"),
        properties=dict(value="20",type="text")
    ),
    dict(
        input='Input', 
        name='n_top_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['n_top_genes'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='subset', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['subset'] or config.show_parameters"),
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
    plot = [
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['plot']['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='show', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['show'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='save', 
        description="typing.Union[bool, str, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['save'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),]
)

def filter_genes_dispersion_f(adata,kwargs):

    sc.pp.filter_genes_dispersion(
        adata,
        flavor=type_formater(kwargs["flavor"],typing.Literal['seurat', 'cell_ranger']),
        min_disp=type_formater(kwargs["min_disp"],typing.Optional[float]),
        max_disp=type_formater(kwargs["max_disp"],typing.Optional[float]),
        min_mean=type_formater(kwargs["min_mean"],typing.Optional[float]),
        max_mean=type_formater(kwargs["max_mean"],typing.Optional[float]),
        n_bins=type_formater(kwargs["n_bins"],int),
        n_top_genes=type_formater(kwargs["n_top_genes"],typing.Optional[int]),
        log=type_formater(kwargs["log"],bool),
        subset=type_formater(kwargs["subset"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def filter_genes_dispersion_plot():

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

config.methods["filter_genes_dispersion"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_genes_dispersion_args,

    function = filter_genes_dispersion_f,

    plot = filter_genes_dispersion_plot,

)
