import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from ..functions import *

from app import app

from .. import config

def args_qc():

    options = node_names(exclude_downstream_from_node=config.selected) 
    optionsthreshold = ["total_counts","n_genes_by_counts"]+[i for i in config.adata.var.columns.values if config.adata.var.dtypes[i] in [bool]]
    
    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options":options
        },
        {
            "input":"MeasureTable",
            "name":"measure",
            "description":"Metrics of quality control to compute",
            "value":"[]",
            "options":optionsthreshold,
            "summary":True
        },
    ]

def f_qc(name_analysis, kwargs):
        
    for l in kwargs["measure"]:
        if "total_counts" == l:
            config.adata.obs[name_analysis+"_total_counts"] = np.array(config.adata.X.sum(axis=1)).reshape(-1)
        elif "n_genes_by_counts" == l:
            config.adata.obs[name_analysis+"_n_genes_by_counts"] = np.array((config.adata.X > 0).sum(axis=1)).reshape(-1)
        else:
            config.adata.obs[name_analysis+"_"+l] = np.array(config.adata.X[:,config.adata.var[l].values].sum(axis=1)).reshape(-1)/config.adata.obs["total_counts"].values
            config.adata.obs[name_analysis+"_"+l] = np.nan_to_num(config.adata.obs[l].values, nan=0)

def rm_qc(name_analysis):

    ls = get_node_parameters(config.selected, str2list=True)["measure"]
    for l in ls:
        name = (name_analysis+"_"+l)
        config.adata.obs.drop(name, axis=1, inplace=True)

def rename_qc(name_analysis, name_new_analysis):

    ls = [(name_analysis+"_"+l) for l in get_node_parameters(name_analysis, str2list=True)["measure"]]
    cols = {}    
    for i in config.adata.obs.columns.values:
        if i not in ls:
            cols[i] = i
        else:
            name = i.split(name_analysis)[-1]
            name = (name_new_analysis+"_"+name)
            cols[i] = name

    config.adata.obs.rename(columns=cols, inplace=True)

def plot_qc(name_analysis):

    l = []
    # if "measure" in config.adata.uns["__interactive__"][name_analysis]["params"].keys():
    #     dic = config.adata.uns["__interactive__"][name_analysis]["params"]["measure"] 
    #     dic = dic[1:-1].split(",")
    #     dic = [i for i in dic if i != '']
    #     for var_selected_data in dic:

    #         # var_selected_data = metric["Measure"]

    #         # hist_values, hist_bins = np.histogram(config.adata.obs[var_selected_data].values, bins=30)
    #         # tallest_bin_height = np.max(hist_values)
    #         # var1_vertical_line_min = go.Scatter(
    #         #     x=[metric["Min"], metric["Min"]],
    #         #     y=[0, tallest_bin_height],
    #         #     mode='lines',
    #         #     name='Min threshold',
    #         #     line=dict(color='red', dash='dash')
    #         # )
    #         # var1_vertical_line_max = go.Scatter(
    #         #     x=[metric["Max"], metric["Max"]],
    #         #     y=[0, tallest_bin_height],
    #         #     mode='lines',
    #         #     name='Max threshold',
    #         #     line=dict(color='green', dash='dash')
    #         # )

    #         l += [
    #                 dbc.Row(
    #                     [dcc.Graph(id="Histogram",
    #                             figure={
    #                                     "data":[
    #                                         go.Histogram(
    #                                             x=config.adata.obs[var_selected_data].values,
    #                                             nbinsx=30,
    #                                             name='Histogram',
    #                                             marker=dict(color='blue'),
    #                                             opacity=0.7
    #                                         ),
    #                                         # var1_vertical_line_min,
    #                                         # var1_vertical_line_max
    #                                     ],
    #                                     "layout":{
    #                                             'title': f'Histogram of {var_selected_data}',
    #                                             'xaxis': {'title': var_selected_data},
    #                                             'yaxis': {'title': 'Count'},
    #                                             'barmode': 'overlay',
    #                                             'width':1500,
    #                                             'height':400,
    #                                     }
    #                                 }
    #                     )],
    #                     justify="center",
    #                     style={'width': '90%', 'margin': 'auto'}
    #                 )
    #             ]
        
    return l