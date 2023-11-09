
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

log1p_args = [ARGINPUT,]

def log1p_f(adata,kwargs):

    sc.pp.log1p(
        adata,
    )
        
    return

config.methods["log1p"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = log1p_args,

    function = log1p_f,

)
