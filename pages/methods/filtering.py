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

def f_filtering(name_analysis, kwargs):
        
    pos = get_node_pos(name_analysis)
    keep = config.graph[pos]["data"]["retained"]
    config.adata = config.adata[keep,:]

    return

def rm_filtering(name_analysis):
        
    return

def rename_filtering(name_analysis):
        
    return

def plot_filtering(name_analysis):

    name = config.active_node_parameters["input"]

    if name == None:
        return []

    ancestors = get_ancestors(name)
    ancestors.append(get_node(name))

    pos = get_node_pos(name_analysis)

    for ancestor in ancestors:
        if "filter" in ancestor["data"]:
            if type(ancestor["data"]["filter"]) == dict:
                for i,j in ancestor["data"]["filter"].items():
                    total_removed = np.zeros_like(j)
                    break
            else:
                j = ancestor["data"]["filter"]
                total_removed = np.zeros_like(j)
                break

            config.graph[pos]["data"]["retained"] = np.array(total_removed)==False

    removed = []
    for ancestor in ancestors:
        if "filter" in ancestor["data"]:
            if type(ancestor["data"]["filter"]) == dict:
                for i,j in ancestor["data"]["filter"].items():
                    removed.append({"name":name_analysis+" "+i,"removed":sum(j)*100/len(j),"retained":sum(np.array(j)==False)*100/len(j)})
                    total_removed += j
            else:
                j = ancestor["data"]["filter"]
                removed.append({"name":ancestor["data"]["id"],"removed":sum(j)*100/len(j),"retained":sum(np.array(j)==False)*100/len(j)})
                total_removed += j

    if len(removed) != 0:
        total_removed = total_removed > 0
        removed.append({"name":"TOTAL","removed":sum(total_removed)*100/len(j),"retained":sum(np.array(total_removed)==False)*100/len(j)})

        config.graph[pos]["data"]["retained"] = np.array(total_removed)==False

    return [
        dash_table.DataTable(
            columns=[
                {"name":"name","id":"name"},
                {"name":"removed (%)","id":"removed"},
                {"name":"retained (%)","id":"retained"}
            ],
            data=removed
        )
    ]

@app.callback(
    dash.Output("analysis_plot","children",allow_duplicate=True),
    dash.Input("analysis_input","value"),
    prevent_initial_call=True
)
def update(_):

    prevent_race("filtering",computed=False)

    return plot_filtering(config.selected)

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