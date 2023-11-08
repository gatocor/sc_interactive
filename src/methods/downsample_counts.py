
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

downsample_counts_args = [ARGINPUT,]

def downsample_counts_f(adata,kwargs):

    sc.pp.downsample_counts(
        adata,
    )
        
    return

config.methods["downsample_counts"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = downsample_counts_args,

    function = downsample_counts_f,

)
