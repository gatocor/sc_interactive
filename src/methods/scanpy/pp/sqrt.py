
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

sqrt_args = [ARGINPUT,
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='chunked', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['chunked'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='chunk_size', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['chunk_size'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),]

def sqrt_f(adata,kwargs):

    sc.pp.sqrt(
        config.adata,
        copy=type_formater(kwargs["copy"],bool),
        chunked=type_formater(kwargs["chunked"],bool),
        chunk_size=type_formater(kwargs["chunk_size"],typing.Optional[int]),
    )
        
    return

config.methods["sqrt"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = sqrt_args,

    function = sqrt_f,

    docs = sc.pp.sqrt.__doc__

)
