
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

recipe_zheng17_args = [ARGINPUT,
    dict(
        input='Input', 
        name='n_top_genes', 
        description="<class 'int'>", 
        visible=dict(function="str(1000)!=config.active_node_parameters['n_top_genes'] or config.show_parameters"),
        properties=dict(value="1000",type="text")
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
        name='plot', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['plot'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),]

def recipe_zheng17_f(adata,kwargs):

    sc.pp.recipe_zheng17(
        adata,
        n_top_genes=type_formater(kwargs["n_top_genes"],int),
        log=type_formater(kwargs["log"],bool),
        plot=type_formater(kwargs["plot"],bool),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["recipe_zheng17"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = recipe_zheng17_args,

    function = recipe_zheng17_f,

    docs = sc.pp.recipe_zheng17.__doc__

)
