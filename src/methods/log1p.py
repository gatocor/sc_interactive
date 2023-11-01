import scanpy as sc

from general import *

log1p_args = {

    "execution" : [
        ARGINPUT,
    ],

    "postexecution" : [],

    "plot" : [],

}

def log1p_f(adata, kwargs):
        
    sc.pp.log1p(config.adata, base=np.e)

def log1p_plot():

    return []

config.methods["log1p"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":True,
    },

    "args": log1p_args,

    "function": log1p_f,

    "plot": log1p_plot,

}