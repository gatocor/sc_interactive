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

def leiden_args(adata):

    options = []
    if "__interactive__" in adata.uns.keys():
        options = [i for i,j in adata.uns["__interactive__"].items() if j["type"]=="Neighbors"] 
    
    return [
        "Leiden",
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
        {
            "input":"Input",
            "name":"n_iterations",
            "description":"How many iterations of the Leiden clustering algorithm to perform. Positive values above 2 define the total number of iterations to perform, -1 has the algorithm run until it reaches its optimal clustering.",
            "value":-1,
            "type":"number"
        },
        # partition_type=None
    ]

def f_leiden(adata, name_analysis, **kwargs):

    sc.tl.leiden(
        adata,
        neighbors_key=kwargs["input"],
        resolution=kwargs["resolution"],
        random_state=kwargs["random_state"], 
        directed=kwargs["directed"], 
        use_weights=kwargs["use_weights"], 
        n_iterations=kwargs["n_iterations"], 
    )

def make_leiden_plots1(adata, name_analysis):

    return []

def make_leiden_plots2(adata, name_analysis):

    return []
