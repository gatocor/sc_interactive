
import numpy
from numpy import inf
import scanpy as sc
import louvain
import scipy
import leidenalg

from general import *

paga_expression_entropies_args = [ARGINPUT,]

def paga_expression_entropies_f(adata,kwargs):

    sc.tl.paga_expression_entropies(
        config.adata,
    )
        
    return

config.methods["paga_expression_entropies"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = paga_expression_entropies_args,

    function = paga_expression_entropies_f,

    docs = sc.tl.paga_expression_entropies.__doc__

)
