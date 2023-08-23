import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash
import inspect

from app import app

from .. import config

def louvain_args(adata):

    options = []
    if "__interactive__" in adata.uns.keys():
        options = [i for i,j in adata.uns["__interactive__"].items() if j["type"]=="Neighbors"] 
    
    return [
        "louvain",
            {
            "input":"Dropdown",
            "name":"input",
            "description":"Neighbor representation to use for clustering.",
            "value":"Raw",
            "clearable":False,
            "options":options
        },
        {
            "input":"Input",
            "name":"resolution",
            "description":"A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. Set to None if overriding partition_type to one that doesn’t accept a resolution_parameter.",
            "value":1,
            "type":"number"
        },
        # restrict_to=None, 
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change the initialization of the optimization.",
            "value":0,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"flavor",
            "description":"Choose between to packages for computing the clustering. 'vtraag' is much more powerful, and the default.",
            "value":'vtraag',
            "clearable":False,
            "options":['vtraag', 'igraph', 'rapids']
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

def f_louvain(adata, name_analysis, **kwargs):

    sc.tl.louvain(
        adata,
        neighbors_key=kwargs["input"],
        resolution=kwargs["resolution"],
        random_state=kwargs["random_state"], 
        flavor=kwargs["flavor"], 
        directed=kwargs["directed"], 
        use_weights=kwargs["use_weights"], 
    )

def make_louvain_plots1(adata, name_analysis):

    return []

def make_louvain_plots2(adata, name_analysis):

    return []
