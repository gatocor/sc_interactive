import scanpy as sc
import pandas as pd
import os

# Initialize the global variable
analysis_folder = None
h5ad_file = "Raw.h5ad"
adata = None
graph = []
report = None
show_parameters = False
show_plot = False

active_node_parameters = {}
active_plot_parameters = {}
max_x = 0
selected = None
selected_plot = None
figure = None
add_method = None

methods = {
    "Raw" : {
        "properties" : {
            "type" : "Raw",
        },

        "args": [],

        "plots":[],

        "docs":[]
    }
}
methods_plot = {}