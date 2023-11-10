
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

diffmap_args = [ARGINPUT,
    dict(
        input='Input', 
        name='n_comps', 
        description="<class 'int'>", 
        visible=dict(function="str(15)!=config.active_node_parameters['n_comps'] or config.show_parameters"),
        properties=dict(value="15",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['neighbors_key'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),]

def diffmap_f(adata,kwargs):

    sc.tl.diffmap(
        config.adata,
        n_comps=type_formater(kwargs["n_comps"],int),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["diffmap"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = diffmap_args,

    function = diffmap_f,

    docs = sc.tl.diffmap.__doc__

)
