
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

filter_rank_genes_groups_args = [ARGINPUT,]

def filter_rank_genes_groups_f(adata,kwargs):

    sc.tl.filter_rank_genes_groups(
        adata,
    )
        
    return

config.methods["filter_rank_genes_groups"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = filter_rank_genes_groups_args,

    function = filter_rank_genes_groups_f,

    docs = sc.tl.filter_rank_genes_groups.__doc__

)
