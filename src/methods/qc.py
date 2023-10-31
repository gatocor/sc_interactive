import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash
from scipy.stats import mode
import plotly.express as px

from general import *

args = {

    "execution" : [
        ARGINPUT,
        ARGMATRIX,
        ARGMATRIX_KEY,
        ARGBATCH,
        {
            "input":"AgTable",
            "name":"measure",
            "description":"Metrics of quality control to compute",
            "properties":{
                "header":[
                    { "headerName": "Name", "field":"name", "editable": True },
                    { "headerName": ".var", "field":"var", "editable": True,
                    "cellEditor": "agSelectCellEditor",
                    "cellEditorParams": {"values": {"function":"[None]+list(config.adata.var.columns.values)"}},
                    },
                    { "headerName": "Pattern", "field":"pattern", "editable": True },
                    { "headerName": "Style", "field":"style", "editable": True,
                    "cellEditor": "agSelectCellEditor",
                    "cellEditorParams": {"values": ['counts','n_expressed_genes','proportion']},
                    }
                ],
                "data":[
                    {
                        "name":"counts", "var":"", "pattern":"", "style":"counts",
                    },
                    {
                        "name":"n_expressed_genes", "var":"", "pattern":"", "style":"n_expressed_genes",
                    }
                ],
            },
            "addRows":{"name":"", "var":"", "pattern":"", "style":"counts"},
            "deleteRows": True                
        },
    ],

    "postexecution" : [
        {
            "input":"AgTable",
            "name":"thresholds",
            "description":"Thresholds to apply to the data",
            "properties":{
                "header":[
                    { "headerName": "Metric", "field":"metric", "editable": False },
                    { "headerName": "Batch", "field":"batch", "editable": False },
                    { "headerName": "Min", "field":"min", "editable": True },
                    { "headerName": "Max", "field":"max", "editable": True },
                    { "headerName": "Genes", "field":"genes", "editable": False },
                ],
                "data":{"function":"qc_data()"},
            }
        },
    ],

    "plot" : [
        {
            "input":"Dropdown",
            "name":"plot_style",
            "description":"Chose between violin or 2D scatterplot.",
            "properties" : {
                "value":"violin",
                "clearable":False,
                "options":["violin","scatter"],
            },
        },
        {
            "input":"Dropdown",
            "name":"plot_x",
            "description":"Genes to use to for the dimensionality reduction.",
            "properties":{
                "value":{"function": "qc_measures()[0]"},
                "clearable":False,
                "options":{"function": "qc_measures()"},
            },
            "visible":{"function": "get_node(config.selected)['data']['plot']['plot_style'] == 'scatter'"}
        },
        {
            "input":"Dropdown",
            "name":"plot_y",
            "description":"Genes to use to for the dimensionality reduction.",
            "properties":{
                "value":{"function": "qc_measures()[0]"},
                "clearable":False,
                "options":{"function": "qc_measures()"},
            },
            # "visible":{"function": "get_node(config.selected)['data']['plot']['plot_style'] == 'scatter'"}
        },
        {
            "input":"Dropdown",
            "name":"plot_color",
            "description":"Genes to use to for the dimensionality reduction.",
            "properties":{
                "value":None,
                "clearable":True,
                "options":{"function": "list(config.adata.obs.columns.values)"},
            },
            "visible":{"function": "get_node(config.selected)['data']['plot']['plot_style'] == 'scatter'"}
        },
    ]

}

def qc_measures():

    selected = config.selected
    node = get_node(selected)["data"]["parameters"]

    return [f"{selected}--{i['name']}" for i in node["measure"]]

def qc_data():

    adata = config.adata

    parameters = get_node(config.selected)["data"]["parameters"]

    measures = [i["name"] for i in parameters["measure"]]

    batch = [""]
    if parameters["batch"] != None:
        batch = np.unique(adata.obs[parameters["batch"]].values)

    data = []
    for i in parameters["measure"]:
        genes = qc_get_genes(i.copy())
        if i["name"] in measures:
            j = f"{config.selected}--{i['name']}"
            for b in batch:
                data.append({"metric":i["name"], "batch":str(b),"min":str(adata.obs[j].min()),"max":str(adata.obs[j].max()),"genes":str(genes)})

    return data

def qc_get_genes(b):
    if b["var"] not in  [''] and b["pattern"] != '':
        v =  [i for i in config.adata.var[b["var"]].values if i.startswith(b["pattern"])]
        if len(v) != config.adata.shape[1]:
            genes = v
        else:
            genes = "all"
    else:
        genes = "all"

    return genes

def qc_f(adata, kwargs):
        
    X = get_matrix(kwargs["matrix"], kwargs["matrix_key"])

    d = {"obs":{}}
    for l in kwargs["measure"]:
        genes = qc_get_genes(l)
        if genes == "all":
            x = range(X.shape[0])
        else:
            x = np.isin(adata.var[l["var"]].values, genes)

        name = f"{config.selected}--{l['name']}"
        if "counts" == l["style"]:
            config.adata.obs[name] = np.array(X[:,x].sum(axis=1)).reshape(-1)
        elif "n_expressed_genes" == l["style"]:
            config.adata.obs[name] = np.array((X[:,x] > 0).sum(axis=1)).reshape(-1)
        elif "proportion" == l["style"]:
            p = np.array(config.adata.X[:,x].sum(axis=1)).reshape(-1)/np.array(config.adata.X.sum(axis=1)).reshape(-1)
            config.adata.obs[name] = p
            
        config.adata.obs[name+"--keep"] = np.ones(X.shape[0], bool)

    return

def qc_reset_lims():

    batch = get_node(config.selected)["data"]["parameters"]["batch"]
    for l in get_node(config.selected)["data"]["parameters"]["thresholds"]:
        name = get_name(l["metric"])
        name_keep = get_name(l["metric"])+"--keep"
        if batch != None:
            sub = config.adata.obs[batch] == l["batch"]
            config.adata.obs[name_keep].values[sub]  = config.adata.obs[name].values[sub] >= float(l["min"])
            config.adata.obs[name_keep].values[sub]  = config.adata.obs[name].values[sub] <= float(l["max"])
        else:
            config.adata.obs[name_keep]  = config.adata.obs[name].values >= float(l["min"])
            config.adata.obs[name_keep]  = config.adata.obs[name].values <= float(l["max"])

def qc_plot():

    #Reset threshold lims
    qc_reset_lims()

    params = get_node(config.selected)["data"]["parameters"]
    plot_params = get_node(config.selected)["data"]["plot"]

    if plot_params["plot_style"] == "violin":

        if params['batch'] != None:
            x = config.adata.obs[params['batch']].values
        else:
            x = [0 for i in range(config.adata.shape[0])]

        y = config.adata.obs[plot_params['plot_y']].values

        thresholds = get_from_table(params['thresholds'],filter={"metric":[plot_params['plot_y'].split("--")[1]]},properties={"batch":str,"min":float,"max":float})
        if "" in thresholds['batch']:
            thresholds['batch'] = [0 for i in thresholds['batch']]

        fig = px.violin(x=x,
                        y=y
              )
        fig.add_traces(list(
            px.scatter(x=thresholds['batch'],
                       y=thresholds['min'],
                       color_discrete_sequence=["red"], 
                       symbol_sequence=["triangle-up-dot"]
            ).select_traces()
        ))
        fig.add_traces(list(
            px.scatter(x=thresholds['batch'],
                       y=thresholds['max'],
                       color_discrete_sequence=["purple"], 
                       symbol_sequence=["triangle-down-dot"]
            ).select_traces()
        ))

        fig.layout["xaxis"]["title"] = params['batch']
        fig.layout["yaxis"]["title"] = plot_params['plot_y']

        return dcc.Graph(figure=fig)
        
    else:

        x = config.adata.obs[plot_params['plot_x']].values
        y = config.adata.obs[plot_params['plot_y']].values

        if plot_params['plot_color'] != None:
            color = config.adata.obs[plot_params['plot_color']].values

            fig = px.scatter(x=x,y=y,color=color)
            fig.layout["legend"]["title"]["text"] = plot_params['plot_color']
        else:
            fig = px.scatter(x=x,y=y)

        fig.layout["xaxis"]["title"] = plot_params['plot_x']
        fig.layout["yaxis"]["title"] = plot_params['plot_y']

        return plot_center(dcc.Graph(figure=fig, style={"width": "60vw", "height": "50vw"}))

config.methods["qc"] = {

    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": args,

    "function":qc_f,

    "plot":qc_plot,

}