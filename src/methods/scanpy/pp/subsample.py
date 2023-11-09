
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

subsample_args = [ARGINPUT,
    dict(
        input='Input', 
        name='fraction', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['fraction'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_obs', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['n_obs'] or config.show_parameters"),
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

def subsample_f(adata,kwargs):

    sc.pp.subsample(
        adata,
        fraction=type_formater(kwargs["fraction"],typing.Optional[float]),
        n_obs=type_formater(kwargs["n_obs"],typing.Optional[int]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["subsample"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = subsample_args,

    function = subsample_f,

    docs = sc.pp.subsample.__doc__

)
