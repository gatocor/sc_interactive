
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

highly_variable_genes_args = [ARGINPUT,
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['layer'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_top_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['n_top_genes'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='min_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(0.5)!=config.active_node_parameters['min_disp'] or config.show_parameters"),
        properties=dict(value="0.5",type="text")
    ),
    dict(
        input='Input', 
        name='max_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(inf)!=config.active_node_parameters['max_disp'] or config.show_parameters"),
        properties=dict(value="inf",type="text")
    ),
    dict(
        input='Input', 
        name='min_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(0.0125)!=config.active_node_parameters['min_mean'] or config.show_parameters"),
        properties=dict(value="0.0125",type="text")
    ),
    dict(
        input='Input', 
        name='max_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(3)!=config.active_node_parameters['max_mean'] or config.show_parameters"),
        properties=dict(value="3",type="text")
    ),
    dict(
        input='Input', 
        name='span', 
        description="typing.Optional[float]", 
        visible=dict(function="str(0.3)!=config.active_node_parameters['span'] or config.show_parameters"),
        properties=dict(value="0.3",type="text")
    ),
    dict(
        input='Input', 
        name='n_bins', 
        description="<class 'int'>", 
        visible=dict(function="str(20)!=config.active_node_parameters['n_bins'] or config.show_parameters"),
        properties=dict(value="20",type="text")
    ),
    dict(
        input='Input', 
        name='flavor', 
        description="typing.Literal['seurat', 'cell_ranger', 'seurat_v3']", 
        visible=dict(function="'seurat'!=eval(config.active_node_parameters['flavor']) or config.show_parameters"),
        properties=dict(value="'seurat'",type="text")
    ),
    dict(
        input='Input', 
        name='subset', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['subset'] or config.show_parameters"),
        properties=dict(value="False",type="text")
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
        name='batch_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['batch_key'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='check_values', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['check_values'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),]

def highly_variable_genes_f(adata,kwargs):

    sc.pp.highly_variable_genes(
        adata,
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        n_top_genes=type_formater(kwargs["n_top_genes"],typing.Optional[int]),
        min_disp=type_formater(kwargs["min_disp"],typing.Optional[float]),
        max_disp=type_formater(kwargs["max_disp"],typing.Optional[float]),
        min_mean=type_formater(kwargs["min_mean"],typing.Optional[float]),
        max_mean=type_formater(kwargs["max_mean"],typing.Optional[float]),
        span=type_formater(kwargs["span"],typing.Optional[float]),
        n_bins=type_formater(kwargs["n_bins"],int),
        flavor=type_formater(kwargs["flavor"],typing.Literal['seurat', 'cell_ranger', 'seurat_v3']),
        subset=type_formater(kwargs["subset"],bool),
        inplace=type_formater(kwargs["inplace"],bool),
        batch_key=type_formater(kwargs["batch_key"],typing.Optional[str]),
        check_values=type_formater(kwargs["check_values"],bool),
    )
        
    return

config.methods["highly_variable_genes"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = highly_variable_genes_args,

    function = highly_variable_genes_f,

)
