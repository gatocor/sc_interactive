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

def args_filtering():

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

def f_filtering(adata, name_analysis, **kwargs):
        
    return

def rm_filtering(name_analysis):
        
    return

def rename_filtering(name_analysis):
        
    return

def plot_filtering(name_analysis):

    return []

# @app.callback(
#     dash.Output("filtering_measure","data"),
#     dash.Input("filtering_measure_button","n_clicks"),
#     dash.State("filtering_measure_dropdown","value"),
#     dash.State("filtering_measure","data"),
#     prevent_initial_callback=True
# )
# def add_metric(_,value,data):

#     values = [i["Measure"] for i in data]
#     if value != None and value not in values:
#         data.append({"Measure":value})

#     return data

# @app.callback(
#     dash.Output("filtering_thresholds","data"),
#     dash.Input("filtering_input","value"),
#     prevent_initial_callback=True
# )
# def show_table(value):

#     d = []
#     if value != None:
#         l = [i for i,j in config.adata.uns["__interactive__"].items() if j["type"]=="QC"]
#         d = []
#         for i in l:
#             for j in config.adata.uns["__interactive__"][i]["params"]["measure"][1:-1].split(","):
#                 d.append({"Measure":j,"Min":0,"Max":0,"#bins":30})

#     if "params" not in config.adata.uns["__interactive__"][config.name_analysis].keys():
#         config.adata.uns["__interactive__"][config.name_analysis]["params"] = {}

#     config.adata.uns["__interactive__"][config.name_analysis]["params"]["thresholds"] = d

#     return d