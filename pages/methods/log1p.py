import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from app import app

from .. import config

def log1p_args(adata):

    options = []
    if "__interactive__" in adata.uns.keys():
        options = ["Raw"]+[i[2:] for i in adata.obsm.keys()] 
    
    return [
        "Log1p",
            {
            "input":"Dropdown",
            "name":"input",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":"Raw",
            "clearable":False,
            "options":options
        },

    ]

def f_log1p(adata, name_analysis, **kwargs):
        
    if kwargs["input"] == "Raw":
        adata.obsm["X_"+name_analysis] = np.log1p(adata.X)
    else:
        adata.obsm["X_"+name_analysis] = np.log1p(adata.obsm["X_"+kwargs["input"]])

def make_log1p_plots1(adata, name_analysis):

    return []

def make_log1p_plots2(adata, name_analysis):

    return []
