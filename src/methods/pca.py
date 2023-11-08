
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

pca_args = [ARGINPUT,
    dict(
        input='Input', 
        name='n_comps', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['n_comps'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='zero_center', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(True)!=config.active_node_parameters['zero_center'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='svd_solver', 
        description="<class 'str'>", 
        visible=dict(function="'arpack'!=eval(config.active_node_parameters['svd_solver']) or config.show_parameters"),
        properties=dict(value="'arpack'",type="text")
    ),
    dict(
        input='Input', 
        name='random_state', 
        description="typing.Union[str, int, numpy.random.mtrand.RandomState]", 
        visible=dict(function="str(0)!=config.active_node_parameters['random_state'] or config.show_parameters"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='return_info', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['return_info'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='use_highly_variable', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_node_parameters['use_highly_variable'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='dtype', 
        description="<class 'str'>", 
        visible=dict(function="'float32'!=eval(config.active_node_parameters['dtype']) or config.show_parameters"),
        properties=dict(value="'float32'",type="text")
    ),
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

def pca_f(adata,kwargs):

    sc.tl.pca(
        adata,
        n_comps=type_formater(kwargs["n_comps"],typing.Optional[int]),
        zero_center=type_formater(kwargs["zero_center"],typing.Optional[bool]),
        svd_solver=type_formater(kwargs["svd_solver"],str),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        return_info=type_formater(kwargs["return_info"],bool),
        use_highly_variable=type_formater(kwargs["use_highly_variable"],typing.Optional[bool]),
        dtype=type_formater(kwargs["dtype"],str),
        copy=type_formater(kwargs["copy"],bool),
        chunked=type_formater(kwargs["chunked"],bool),
        chunk_size=type_formater(kwargs["chunk_size"],typing.Optional[int]),
    )
        
    return

config.methods["pca"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = pca_args,

    function = pca_f,

)
