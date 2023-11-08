
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

combat_args = [ARGINPUT,
    dict(
        input='Input', 
        name='key', 
        description="<class 'str'>", 
        visible=dict(function="'batch'!=eval(config.active_node_parameters['key']) or config.show_parameters"),
        properties=dict(value="'batch'",type="text")
    ),
    dict(
        input='Input', 
        name='covariates', 
        description="typing.Optional[typing.Collection[str]]", 
        visible=dict(function="str(None)!=config.active_node_parameters['covariates'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='inplace', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['inplace'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),]

def combat_f(adata,kwargs):

    sc.pp.combat(
        adata,
        key=type_formater(kwargs["key"],str),
        covariates=type_formater(kwargs["covariates"],typing.Optional[typing.Collection[str]]),
        inplace=type_formater(kwargs["inplace"],bool),
    )
        
    return

config.methods["combat"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = combat_args,

    function = combat_f,

)
