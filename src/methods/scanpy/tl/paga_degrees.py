
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

paga_degrees_args = [ARGINPUT,]

def paga_degrees_f(adata,kwargs):

    sc.tl.paga_degrees(
        adata,
    )
        
    return

config.methods["paga_degrees"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = paga_degrees_args,

    function = paga_degrees_f,

    docs = sc.tl.paga_degrees.__doc__

)
