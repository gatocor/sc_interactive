
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

umap_args = [ARGINPUT,
    dict(
        input='Input', 
        name='min_dist', 
        description="<class 'float'>", 
        visible=dict(function="str(0.5)!=config.active_node_parameters['min_dist'] or config.show_parameters"),
        properties=dict(value="0.5",type="text")
    ),
    dict(
        input='Input', 
        name='spread', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['spread'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='n_components', 
        description="<class 'int'>", 
        visible=dict(function="str(2)!=config.active_node_parameters['n_components'] or config.show_parameters"),
        properties=dict(value="2",type="text")
    ),
    dict(
        input='Input', 
        name='maxiter', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['maxiter'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='alpha', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['alpha'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='gamma', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['gamma'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='negative_sample_rate', 
        description="<class 'int'>", 
        visible=dict(function="str(5)!=config.active_node_parameters['negative_sample_rate'] or config.show_parameters"),
        properties=dict(value="5",type="text")
    ),
    dict(
        input='Input', 
        name='init_pos', 
        description="typing.Union[typing.Literal['paga', 'spectral', 'random'], numpy.ndarray, str]", 
        visible=dict(function="'spectral'!=eval(config.active_node_parameters['init_pos']) or config.show_parameters"),
        properties=dict(value="'spectral'",type="text")
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
        name='a', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['a'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='b', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['b'] or config.show_parameters"),
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
        name='method', 
        description="typing.Literal['umap', 'rapids']", 
        visible=dict(function="'umap'!=eval(config.active_node_parameters['method']) or config.show_parameters"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['neighbors_key'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),]

def umap_f(adata,kwargs):

    sc.tl.umap(
        adata,
        min_dist=type_formater(kwargs["min_dist"],float),
        spread=type_formater(kwargs["spread"],float),
        n_components=type_formater(kwargs["n_components"],int),
        maxiter=type_formater(kwargs["maxiter"],typing.Optional[int]),
        alpha=type_formater(kwargs["alpha"],float),
        gamma=type_formater(kwargs["gamma"],float),
        negative_sample_rate=type_formater(kwargs["negative_sample_rate"],int),
        init_pos=type_formater(kwargs["init_pos"],typing.Union[typing.Literal['paga', 'spectral', 'random'], numpy.ndarray, str]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        a=type_formater(kwargs["a"],typing.Optional[float]),
        b=type_formater(kwargs["b"],typing.Optional[float]),
        copy=type_formater(kwargs["copy"],bool),
        method=type_formater(kwargs["method"],typing.Literal['umap', 'rapids']),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
    )
        
    return

config.methods["umap"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = umap_args,

    function = umap_f,

)
