
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

tsne_args = [ARGINPUT,
    dict(
        input='Input', 
        name='n_pcs', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['n_pcs'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_rep', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['use_rep'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='perplexity', 
        description="typing.Union[float, int]", 
        visible=dict(function="str(30)!=config.active_node_parameters['perplexity'] or config.show_parameters"),
        properties=dict(value="30",type="text")
    ),
    dict(
        input='Input', 
        name='early_exaggeration', 
        description="typing.Union[float, int]", 
        visible=dict(function="str(12)!=config.active_node_parameters['early_exaggeration'] or config.show_parameters"),
        properties=dict(value="12",type="text")
    ),
    dict(
        input='Input', 
        name='learning_rate', 
        description="typing.Union[float, int]", 
        visible=dict(function="str(1000)!=config.active_node_parameters['learning_rate'] or config.show_parameters"),
        properties=dict(value="1000",type="text")
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
        name='use_fast_tsne', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['use_fast_tsne'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='n_jobs', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['n_jobs'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='metric', 
        description="<class 'str'>", 
        visible=dict(function="'euclidean'!=eval(config.active_node_parameters['metric']) or config.show_parameters"),
        properties=dict(value="'euclidean'",type="text")
    ),]

def tsne_f(adata,kwargs):

    sc.tl.tsne(
        adata,
        n_pcs=type_formater(kwargs["n_pcs"],typing.Optional[int]),
        use_rep=type_formater(kwargs["use_rep"],typing.Optional[str]),
        perplexity=type_formater(kwargs["perplexity"],typing.Union[float, int]),
        early_exaggeration=type_formater(kwargs["early_exaggeration"],typing.Union[float, int]),
        learning_rate=type_formater(kwargs["learning_rate"],typing.Union[float, int]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        use_fast_tsne=type_formater(kwargs["use_fast_tsne"],bool),
        n_jobs=type_formater(kwargs["n_jobs"],typing.Optional[int]),
        copy=type_formater(kwargs["copy"],bool),
        metric=type_formater(kwargs["metric"],str),
    )
        
    return

config.methods["tsne"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = tsne_args,

    function = tsne_f,

    docs = sc.tl.tsne.__doc__

)
