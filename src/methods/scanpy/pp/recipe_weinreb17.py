
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

recipe_weinreb17_args = [ARGINPUT,
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['log'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='mean_threshold', 
        description="<class 'float'>", 
        visible=dict(function="str(0.01)!=config.active_node_parameters['mean_threshold'] or config.show_parameters"),
        properties=dict(value="0.01",type="text")
    ),
    dict(
        input='Input', 
        name='cv_threshold', 
        description="<class 'int'>", 
        visible=dict(function="str(2)!=config.active_node_parameters['cv_threshold'] or config.show_parameters"),
        properties=dict(value="2",type="text")
    ),
    dict(
        input='Input', 
        name='n_pcs', 
        description="<class 'int'>", 
        visible=dict(function="str(50)!=config.active_node_parameters['n_pcs'] or config.show_parameters"),
        properties=dict(value="50",type="text")
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

def recipe_weinreb17_f(adata,kwargs):

    sc.pp.recipe_weinreb17(
        adata,
        log=type_formater(kwargs["log"],bool),
        mean_threshold=type_formater(kwargs["mean_threshold"],float),
        cv_threshold=type_formater(kwargs["cv_threshold"],int),
        n_pcs=type_formater(kwargs["n_pcs"],int),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["recipe_weinreb17"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = recipe_weinreb17_args,

    function = recipe_weinreb17_f,

)
