
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

embedding_density_args = [ARGINPUT,
    dict(
        input='Input', 
        name='basis', 
        description="<class 'str'>", 
        visible=dict(function="'umap'!=eval(config.active_node_parameters['basis']) or config.show_parameters"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='groupby', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['groupby'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='components', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_node_parameters['components'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),]

def embedding_density_f(adata,kwargs):

    sc.tl.embedding_density(
        adata,
        basis=type_formater(kwargs["basis"],str),
        groupby=type_formater(kwargs["groupby"],typing.Optional[str]),
        key_added=type_formater(kwargs["key_added"],typing.Optional[str]),
        components=type_formater(kwargs["components"],typing.Union[str, typing.Sequence[str]]),
    )
        
    return

config.methods["embedding_density"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = embedding_density_args,

    function = embedding_density_f,

    docs = sc.tl.embedding_density.__doc__

)
