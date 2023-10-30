import scanpy as sc
import pandas as pd
import os

# Initialize the global variable
analysis_folder = None
h5ad_file = "Raw.h5ad"
adata = None
graph = []
report = None

active_node_parameters = {}
max_x = 0
selected = 'Raw'

methods = {
    "Raw" : {
        "properties" : {
            "type" : "Raw",
            "make_new_h5ad" : True
        },

        "args": {
            "execution" : [],
            "postexecution" : [],
            "plot" : []
        } 
    }
}