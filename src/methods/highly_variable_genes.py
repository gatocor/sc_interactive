
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

highly_variable_genes_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['layer'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='min_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(0.5)!=get_node(config.selected)['data']['parameters']['min_disp'] or config.show_parameters"),
        properties=dict(value="0.5",type="text")
    ),
    dict(
        input='Input', 
        name='max_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(inf)!=get_node(config.selected)['data']['parameters']['max_disp'] or config.show_parameters"),
        properties=dict(value="inf",type="text")
    ),
    dict(
        input='Input', 
        name='min_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(0.0125)!=get_node(config.selected)['data']['parameters']['min_mean'] or config.show_parameters"),
        properties=dict(value="0.0125",type="text")
    ),
    dict(
        input='Input', 
        name='max_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(3)!=get_node(config.selected)['data']['parameters']['max_mean'] or config.show_parameters"),
        properties=dict(value="3",type="text")
    ),
    dict(
        input='Input', 
        name='span', 
        description="typing.Optional[float]", 
        visible=dict(function="str(0.3)!=get_node(config.selected)['data']['parameters']['span'] or config.show_parameters"),
        properties=dict(value="0.3",type="text")
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
        name='flavor', 
        description="typing.Literal['seurat', 'cell_ranger', 'seurat_v3']", 
        visible=dict(function="'seurat'!=eval(get_node(config.selected)['data']['parameters']['flavor']) or config.show_parameters"),
        properties=dict(value="'seurat'",type="text")
    ),
    dict(
        input='Input', 
        name='subset', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['subset'] or config.show_parameters"),
        properties=dict(value="False",type="text")
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
        name='batch_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['batch_key'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='check_values', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['check_values'] or config.show_parameters"),
        properties=dict(value="True",type="text")
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
    ),
    dict(
        input='Input', 
        name='highly_variable_genes', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['plot']['highly_variable_genes'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),]
)

def highly_variable_genes_f(adata,kwargs):

    sc.pp.highly_variable_genes(
        adata,
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        n_top_genes=type_formater(kwargs["n_top_genes"],typing.Optional[int]),
        min_disp=type_formater(kwargs["min_disp"],typing.Optional[float]),
        max_disp=type_formater(kwargs["max_disp"],typing.Optional[float]),
        min_mean=type_formater(kwargs["min_mean"],typing.Optional[float]),
        max_mean=type_formater(kwargs["max_mean"],typing.Optional[float]),
        span=type_formater(kwargs["span"],typing.Optional[float]),
        n_bins=type_formater(kwargs["n_bins"],int),
        flavor=type_formater(kwargs["flavor"],typing.Literal['seurat', 'cell_ranger', 'seurat_v3']),
        subset=type_formater(kwargs["subset"],bool),
        inplace=type_formater(kwargs["inplace"],bool),
        batch_key=type_formater(kwargs["batch_key"],typing.Optional[str]),
        check_values=type_formater(kwargs["check_values"],bool),
    )
        
    return

def highly_variable_genes_plot():

    kwargs = get_node(config.selected)['data']['plot']
    
    fig = sc.pl.highly_variable_genes(
        config.adata,
        log=type_formater(kwargs["log"],bool),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
        highly_variable_genes=type_formater(kwargs["highly_variable_genes"],bool),
    )


    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["highly_variable_genes"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = highly_variable_genes_args,

    function = highly_variable_genes_f,

    plot = highly_variable_genes_plot,

)
