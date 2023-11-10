
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

scale_args = [ARGINPUT,]

def scale_f(adata,kwargs):

    sc.pp.scale(
        config.adata,
    )
        
    return

config.methods["scale"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = scale_args,

    function = scale_f,

    docs = sc.pp.scale.__doc__

)
