
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

marker_gene_overlap_args = [ARGINPUT,
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
    ),]

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

config.methods["marker_gene_overlap"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = marker_gene_overlap_args,

    function = marker_gene_overlap_f,

)
