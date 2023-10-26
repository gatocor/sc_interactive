import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
from dash import dcc, dash_table
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from ..functions import *
from ..plots import *

from app import app

from .. import config

def args_normalize():

    options = node_names(exclude_downstream_from_node=config.selected) 
    
    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":None,
            "clearable":False,
            "options":options
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

def f_normalize(name_analysis, kwargs):
        
    pos = get_node_pos(name_analysis)

    config.graph[pos]['data']['total_counts_before'] = np.array(config.adata.X.sum(axis=1)).reshape(-1)

    sc.pp.normalize_total(config.adata,
        target_sum=kwargs["target_sum"],
        exclude_highly_expressed=kwargs["exclude_highly_expressed"],
        max_fraction=kwargs["max_fraction"]
    )

    config.graph[pos]['data']['total_counts_after'] = np.array(config.adata.X.sum(axis=1)).reshape(-1)

def rm_normalize(name_analysis):
        
    return

def rename_normalize(name_analysis):
        
    return

def plot_normalize(name_analysis):

    node = get_node(config.selected)
    if node['data']['computed']:

        pos = get_node_pos(name_analysis)

        return [
                dcc.Graph(id="Histogram",
                    figure={
                            "data":[
                                go.Histogram(
                                    x=config.graph[pos]['data']['total_counts_before'],
                                    nbinsx=100,
                                    name='Histogram',
                                    marker=dict(color='blue'),
                                    opacity=0.7
                                ),
                                go.Histogram(
                                    x=config.graph[pos]['data']['total_counts_after'],
                                    nbinsx=100,
                                    name='Histogram',
                                    marker=dict(color='blue'),
                                    opacity=0.7
                                ),
                            ]
                    }
                )
        ]
    
    else:

        return []
