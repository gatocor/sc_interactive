import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash
import scrublet
from scipy.stats import mode

from ..functions import *

from app import app

from .. import config

def args_louvain():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_batch = get_batch_keys()
    
    return [
            {
            "input":"Dropdown",
            "name":"input",
            "description":"Neighbor representation to use for clustering.",
            "value":None,
            "clearable":False,
            "options":options
        },
        {
            "input":"Dropdown",
            "name":"batch_key",
            "description":"str, optional (default: None) Batch key to use. The highly varying will be computed independently in each set of cells separated by batch. If None, use the full dataset.",
            "value":None,
            "clearable":True,
            "options":options_batch,
            "summary":True
        },
        {
            "input":"Input",
            "name":"resolution",
            "description":"A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. Set to None if overriding partition_type to one that doesn’t accept a resolution_parameter.",
            "value":1,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"flavor",
            "description":"Choose between to packages for computing the clustering.",
            "value":"vtraag",
            "clearable":False,
            "options":["vtraag","igraph","rapids"],
            "summary":True
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change the initialization of the optimization.",
            "value":0,
            "type":"number"
        },
        {
            "input":"BooleanSwitch",
            "name":"directed",
            "description":"Whether to treat the graph as directed or undirected.",
            "value":True,
        },
        {
            "input":"BooleanSwitch",
            "name":"use_weights",
            "description":"If True, edge weights from the graph are used in the computation (placing more emphasis on stronger edges).",
            "value":True,
        },
        # partition_type=None
    ]

def f_louvain(name_analysis, kwargs, sub_name, sub):

    node = get_node(kwargs["input"])
    x = config.adata.uns[kwargs["input"]+"_"+str(sub_name)]["connectivities"]
    adata_copy = sc.AnnData(X = np.zeros([x.shape[0],2]))
    adata_copy.uns["neighbors"] = {"connectivities_key":"connectivities","params":{"method":node["data"]["parameters"]["method"]}}
    adata_copy.obsp["connectivities"] = x

    sc.tl.louvain(
        adata_copy,
        neighbors_key=kwargs["input"],
        resolution=kwargs["resolution"],
        random_state=kwargs["random_state"], 
        directed=kwargs["directed"], 
        use_weights=kwargs["use_weights"], 
        flavor=kwargs["flavor"], 
    )

    pos = get_node_pos(name_analysis)
    if "obs" not in config.graph[pos]["data"]:
        config.graph[pos]["data"]["obs"] = {"louvain":np.zeros([config.adata.shape[0]])}
    elif type(config.graph[pos]["data"]["obs"]["louvain"]) == list:
        config.graph[pos]["data"]["obs"] = {"louvain":np.zeros([config.adata.shape[0]])}

    config.graph[pos]["data"]["obs"]["louvain"][sub] = adata_copy.obs["louvain"].values

def rm_louvain(name_analysis):

    return

def rename_louvain(name_analysis, name_new_analysis):

    return

def plot_louvain(name_analysis):

    return []

