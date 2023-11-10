import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

filter_by_obs_args = [
    ARGINPUT,
    dict(
        input='Input', 
        name='obs_key', 
        description="str", 
        visible=True,
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='obs_values', 
        description="typing.Union[str,typing.Collection[str]]", 
        visible=True,
        properties=dict(value="None",type="text")
    ),
]

def filter_by_obs_f(adata,kwargs):
    """
    Filter cells by obs.
    """

    key = type_formater(kwargs["obs_key"],str)
    vals = type_formater(kwargs["obs_values"],typing.Union[str,typing.Collection[str]])

    keep = [i in vals for i in adata.obs[key].values]
    config.adata = config.adata[keep,:]

    return

config.methods["filter_by_obs"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_by_obs_args,

    function = filter_by_obs_f,

    docs = filter_by_obs_f.__doc__

)
