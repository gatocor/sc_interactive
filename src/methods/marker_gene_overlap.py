
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

marker_gene_overlap_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='key', 
        description="<class 'str'>", 
        visible=dict(function="'rank_genes_groups'!=eval(config.active_node_parameters['key']) or config.show_parameters"),
        properties=dict(value="'rank_genes_groups'",type="text")
    ),
    dict(
        input='Input', 
        name='method', 
        description="typing.Literal['overlap_count', 'overlap_coef', 'jaccard']", 
        visible=dict(function="'overlap_count'!=eval(config.active_node_parameters['method']) or config.show_parameters"),
        properties=dict(value="'overlap_count'",type="text")
    ),
    dict(
        input='Input', 
        name='normalize', 
        description="typing.Optional[typing.Literal['reference', 'data']]", 
        visible=dict(function="str(None)!=config.active_node_parameters['normalize'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='top_n_markers', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['top_n_markers'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='adj_pval_threshold', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['adj_pval_threshold'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='key_added', 
        description="<class 'str'>", 
        visible=dict(function="'marker_gene_overlap'!=eval(config.active_node_parameters['key_added']) or config.show_parameters"),
        properties=dict(value="'marker_gene_overlap'",type="text")
    ),
    dict(
        input='Input', 
        name='inplace', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['inplace'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),],
    postexecution = [],
    plot = []
)

def marker_gene_overlap_f(adata,kwargs):

    sc.tl.marker_gene_overlap(
        adata,
        key=type_formater(kwargs["key"],str),
        method=type_formater(kwargs["method"],typing.Literal['overlap_count', 'overlap_coef', 'jaccard']),
        normalize=type_formater(kwargs["normalize"],typing.Optional[typing.Literal['reference', 'data']]),
        top_n_markers=type_formater(kwargs["top_n_markers"],typing.Optional[int]),
        adj_pval_threshold=type_formater(kwargs["adj_pval_threshold"],typing.Optional[float]),
        key_added=type_formater(kwargs["key_added"],str),
        inplace=type_formater(kwargs["inplace"],bool),
    )
        
    return

def marker_gene_overlap_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["marker_gene_overlap"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = marker_gene_overlap_args,

    function = marker_gene_overlap_f,

    plot = marker_gene_overlap_plot,

)
