
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

filter_genes_dispersion_args = [ARGINPUT,
    dict(
        input='Input', 
        name='flavor', 
        description="typing.Literal['seurat', 'cell_ranger']", 
        visible=dict(function="'seurat'!=eval(config.active_node_parameters['flavor']) or config.show_parameters"),
        properties=dict(value="'seurat'",type="text")
    ),
    dict(
        input='Input', 
        name='min_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['min_disp'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_disp', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['max_disp'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='min_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['min_mean'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_mean', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['max_mean'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='n_top_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['n_top_genes'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['log'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='subset', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['subset'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),]

def filter_genes_dispersion_f(adata,kwargs):

    sc.pp.filter_genes_dispersion(
        config.adata,
        flavor=type_formater(kwargs["flavor"],typing.Literal['seurat', 'cell_ranger']),
        min_disp=type_formater(kwargs["min_disp"],typing.Optional[float]),
        max_disp=type_formater(kwargs["max_disp"],typing.Optional[float]),
        min_mean=type_formater(kwargs["min_mean"],typing.Optional[float]),
        max_mean=type_formater(kwargs["max_mean"],typing.Optional[float]),
        n_bins=type_formater(kwargs["n_bins"],int),
        n_top_genes=type_formater(kwargs["n_top_genes"],typing.Optional[int]),
        log=type_formater(kwargs["log"],bool),
        subset=type_formater(kwargs["subset"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["filter_genes_dispersion"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_genes_dispersion_args,

    function = filter_genes_dispersion_f,

    docs = sc.pp.filter_genes_dispersion.__doc__

)
