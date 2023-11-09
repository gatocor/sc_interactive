
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

neighbors_args = [ARGINPUT,
    dict(
        input='Input', 
        name='n_neighbors', 
        description="<class 'int'>", 
        visible=dict(function="str(15)!=config.active_node_parameters['n_neighbors'] or config.show_parameters"),
        properties=dict(value="15",type="text")
    ),
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
        name='knn', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_node_parameters['knn'] or config.show_parameters"),
        properties=dict(value="True",type="text")
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
        name='method', 
        description="typing.Optional[typing.Literal['umap', 'gauss', 'rapids']]", 
        visible=dict(function="'umap'!=eval(config.active_node_parameters['method']) or config.show_parameters"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='metric', 
        description="typing.Union[typing.Literal['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'], typing.Literal['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'], typing.Callable[[numpy.ndarray, numpy.ndarray], float]]", 
        visible=dict(function="'euclidean'!=eval(config.active_node_parameters['metric']) or config.show_parameters"),
        properties=dict(value="'euclidean'",type="text")
    ),
    dict(
        input='Input', 
        name='metric_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=config.active_node_parameters['metric_kwds'] or config.show_parameters"),
        properties=dict(value="{}",type="text")
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
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),]

def neighbors_f(adata,kwargs):

    sc.pp.neighbors(
        adata,
        n_neighbors=type_formater(kwargs["n_neighbors"],int),
        n_pcs=type_formater(kwargs["n_pcs"],typing.Optional[int]),
        use_rep=type_formater(kwargs["use_rep"],typing.Optional[str]),
        knn=type_formater(kwargs["knn"],bool),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        method=type_formater(kwargs["method"],typing.Optional[typing.Literal['umap', 'gauss', 'rapids']]),
        metric=type_formater(kwargs["metric"],typing.Union[typing.Literal['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'], typing.Literal['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'], typing.Callable[[numpy.ndarray, numpy.ndarray], float]]),
        metric_kwds=type_formater(kwargs["metric_kwds"],typing.Mapping[str, typing.Any]),
        key_added=type_formater(kwargs["key_added"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["neighbors"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = neighbors_args,

    function = neighbors_f,

)
