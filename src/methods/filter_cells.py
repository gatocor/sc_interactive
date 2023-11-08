
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

filter_cells_args = [ARGINPUT,
    dict(
        input='Input', 
        name='min_counts', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['min_counts'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='min_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['min_genes'] or config.show_parameters"),
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
        name='max_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['max_genes'] or config.show_parameters"),
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

def filter_cells_f(adata,kwargs):

    sc.pp.filter_cells(
        adata,
        min_counts=type_formater(kwargs["min_counts"],typing.Optional[int]),
        min_genes=type_formater(kwargs["min_genes"],typing.Optional[int]),
        max_counts=type_formater(kwargs["max_counts"],typing.Optional[int]),
        max_genes=type_formater(kwargs["max_genes"],typing.Optional[int]),
        inplace=type_formater(kwargs["inplace"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["filter_cells"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_cells_args,

    function = filter_cells_f,

)
