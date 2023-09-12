import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
from dash import dcc, dash_table
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from ..functions import *

from app import app

from .. import config

def args_log1p():

    options = node_names(exclude_downstream_from_node=config.selected) 
    
    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options":options
        },
    ]

def f_log1p(name_analysis, kwargs):
        
    sc.pp.log1p(config.adata)

    return

def rm_log1p(name_analysis):
        
    return

def rename_log1p(name_analysis):
        
    return

def plot_log1p(name_analysis):

    return []