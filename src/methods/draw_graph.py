
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

draw_graph_args = [ARGINPUT,
    dict(
        input='Input', 
        name='layout', 
        description="typing.Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa']", 
        visible=dict(function="'fa'!=eval(config.active_node_parameters['layout']) or config.show_parameters"),
        properties=dict(value="'fa'",type="text")
    ),
    dict(
        input='Input', 
        name='init_pos', 
        description="typing.Union[str, bool, str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['init_pos'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='root', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['root'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='n_jobs', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['n_jobs'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='adjacency', 
        description="typing.Optional[scipy.sparse._matrix.spmatrix]", 
        visible=dict(function="str(None)!=config.active_node_parameters['adjacency'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='key_added_ext', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['key_added_ext'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='obsp', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['obsp'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),]

def draw_graph_f(adata,kwargs):

    sc.tl.draw_graph(
        adata,
        layout=type_formater(kwargs["layout"],typing.Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa']),
        init_pos=type_formater(kwargs["init_pos"],typing.Union[str, bool, str]),
        root=type_formater(kwargs["root"],typing.Optional[int]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        n_jobs=type_formater(kwargs["n_jobs"],typing.Optional[int]),
        adjacency=type_formater(kwargs["adjacency"],typing.Optional[scipy.sparse._matrix.spmatrix]),
        key_added_ext=type_formater(kwargs["key_added_ext"],typing.Optional[str]),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        obsp=type_formater(kwargs["obsp"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

config.methods["draw_graph"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = draw_graph_args,

    function = draw_graph_f,

)
