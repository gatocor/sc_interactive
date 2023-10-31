import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
from dash import dcc, dash_table
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from general import *

from app import app

args = {
    "execution" : [
            ARGINPUT,
            {
                "input":"AgTable",
                "name":"filtering_thresholds",
                "description":"Metrics of quality control to compute",
                "header":[
                    { "headerName": ".obs", "field":"obs", "editable": True,
                    "cellEditor": "agSelectCellEditor",
                    "cellEditorParams": {"values": {"function":"filtering_names()"}},
                    },
                ],
                "value":{"function":"filtering_dict()"},
                "addRows":{"name":""},
                "deleteRows": True,
                "recomputeAfter": ["input"] 
            },
    ],

    "postexecution" : [],

    "plot" : []
}

def filtering_names():

    if config.active_node_parameters["input"] == None:
        return []
    else:
        return [i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in [bool]]

def filtering_dict():

    if config.active_node_parameters["input"] == None:
        return []
    else:
        return [{"obs":i} for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in [bool] and i.endswith("--keep")]

def filtering_f(adata, inputArgs, kwargs):

    keep = np.ones(adata.X.shape[0], bool)
    for i in kwargs["filtering_thresholds"]:
        name = i["obs"]
        keep = adata.obs[name].values and keep

    d = {
        "keep_cells" : keep
    }

    return d

def plot_filtering():

    l = []
    keep = np.ones(config.adata.X.shape[0], bool)
    for i in config.adata.uns[config.selected]["parameters"]["filtering_thresholds"]:
        name = i["obs"]
        keep = config.adata.obs[name].values and keep
        l.append({"obs":name,"retained":np.mean(config.adata.obs[name].values)})

    l.append({"obs":"TOTAL","retained":np.mean(keep)})

    return [
        ag_table(
            id="filtering_ag_table",
            columnDefs=[
                {"name":"obs","field":"obs"},
                {"name":"retained","field":"retained"}
            ],
            rowData=l
        )
    ]

    name = config.active_node_parameters["input"]

    if name == None:
        return []

    ancestors = get_node_ancestors(name)
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

config.methods["filtering"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":True,
    },

    "args": args,

    "function": filtering_f,

    "plot": plot_filtering,

}