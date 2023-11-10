
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

paga_args = [ARGINPUT,
    dict(
        input='Input', 
        name='groups', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['groups'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_rna_velocity', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['use_rna_velocity'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='model', 
        description="typing.Literal['v1.2', 'v1.0']", 
        visible=dict(function="'v1.2'!=eval(config.active_node_parameters['model']) or config.show_parameters"),
        properties=dict(value="'v1.2'",type="text")
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

def paga_f(adata,kwargs):

    sc.tl.paga(
        config.adata,
        groups=type_formater(kwargs["groups"],typing.Optional[str]),
        use_rna_velocity=type_formater(kwargs["use_rna_velocity"],bool),
        model=type_formater(kwargs["model"],typing.Literal['v1.2', 'v1.0']),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["paga"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = paga_args,

    function = paga_f,

    docs = sc.tl.paga.__doc__

)
