
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

dpt_args = [ARGINPUT,
    dict(
        input='Input', 
        name='n_dcs', 
        description="<class 'int'>", 
        visible=dict(function="str(10)!=config.active_node_parameters['n_dcs'] or config.show_parameters"),
        properties=dict(value="10",type="text")
    ),
    dict(
        input='Input', 
        name='n_branchings', 
        description="<class 'int'>", 
        visible=dict(function="str(0)!=config.active_node_parameters['n_branchings'] or config.show_parameters"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='min_group_size', 
        description="<class 'float'>", 
        visible=dict(function="str(0.01)!=config.active_node_parameters['min_group_size'] or config.show_parameters"),
        properties=dict(value="0.01",type="text")
    ),
    dict(
        input='Input', 
        name='allow_kendall_tau_shift', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['allow_kendall_tau_shift'] or config.show_parameters"),
        properties=dict(value="True",type="text")
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
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),]

def dpt_f(adata,kwargs):

    sc.tl.dpt(
        adata,
        n_dcs=type_formater(kwargs["n_dcs"],int),
        n_branchings=type_formater(kwargs["n_branchings"],int),
        min_group_size=type_formater(kwargs["min_group_size"],float),
        allow_kendall_tau_shift=type_formater(kwargs["allow_kendall_tau_shift"],bool),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["dpt"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = dpt_args,

    function = dpt_f,

)
