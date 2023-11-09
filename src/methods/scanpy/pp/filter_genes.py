
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

filter_genes_args = [ARGINPUT,
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
    ),]

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

config.methods["filter_genes"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_genes_args,

    function = filter_genes_f,

)
