import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from app import app

from .. import config

def normalize_args(adata):

    options = []
    if "__interactive__" in adata.uns.keys():
        options = ["Raw"]+[i[2:] for i in adata.obsm.keys()] 
    
    return [
        "Normalization",
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":"Raw",
            "clearable":False,
            "options":options
        },
        {
            "input":"BooleanSwitch",
            "name":"normalization",
            "description":"If True, normalize the count matrix before PCA analysis.",
            "value":True,
        },
        {
            "input":"Input",
            "name":"target_sum",
            "description":"If None, after normalization, each observation (cell) has a total count equal to the median of total counts for observations (cells) before normalization.",
            "value":None,
            "type":"number"
        },
        {
            "input":"BooleanSwitch",
            "name":"exclude_highly_expressed",
            "description":"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell. A gene is considered highly expressed, if it has more than max_fraction of the total counts in at least one cell. The not-excluded genes will sum up to target_sum.",
            "value":False,
        },
        {
            "input":"Input",
            "name":"max_fraction",
            "description":"If exclude_highly_expressed=True, consider cells as highly expressed that have more counts than max_fraction of the original total counts in at least one cell.",
            "value":0.05,
            "type":"number"
        },
    ]

def f_normalize(adata, name_analysis, **kwargs):
        
    if kwargs["input"] == "Raw":
        adata_copy = sc.AnnData(X=adata.X.copy())
    else:
        adata_copy = sc.AnnData(X=adata.obsm["X_"+kwargs["input"]].copy())

    if kwargs["normalization"]:
        sc.pp.normalize_total(adata_copy,
                              target_sum=kwargs["target_sum"],
                              exclude_highly_expressed=kwargs["exclude_highly_expressed"],
                              max_fraction=kwargs["max_fraction"]
                              )

    adata.obsm["X_"+name_analysis] = adata_copy.X

def make_normalize_plots1(adata, name_analysis):

    return []

def make_normalize_plots2(adata, name_analysis):

    return []
