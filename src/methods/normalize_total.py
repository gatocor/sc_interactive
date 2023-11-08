
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

normalize_total_args = [ARGINPUT,
    dict(
        input='Input', 
        name='target_sum', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['target_sum'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='exclude_highly_expressed', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['exclude_highly_expressed'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='max_fraction', 
        description="<class 'float'>", 
        visible=dict(function="str(0.05)!=config.active_node_parameters['max_fraction'] or config.show_parameters"),
        properties=dict(value="0.05",type="text")
    ),
    dict(
        input='Input', 
        name='key_added', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['key_added'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['layer'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layers', 
        description="typing.Union[typing.Literal['all'], typing.Iterable[str]]", 
        visible=dict(function="str(None)!=config.active_node_parameters['layers'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layer_norm', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['layer_norm'] or config.show_parameters"),
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

def normalize_total_f(adata,kwargs):

    sc.pp.normalize_total(
        adata,
        target_sum=type_formater(kwargs["target_sum"],typing.Optional[float]),
        exclude_highly_expressed=type_formater(kwargs["exclude_highly_expressed"],bool),
        max_fraction=type_formater(kwargs["max_fraction"],float),
        key_added=type_formater(kwargs["key_added"],typing.Optional[str]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        layers=type_formater(kwargs["layers"],typing.Union[typing.Literal['all'], typing.Iterable[str]]),
        layer_norm=type_formater(kwargs["layer_norm"],typing.Optional[str]),
        inplace=type_formater(kwargs["inplace"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["normalize_total"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = normalize_total_args,

    function = normalize_total_f,

)
