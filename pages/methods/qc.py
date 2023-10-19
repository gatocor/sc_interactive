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

from ..arguments import *
from ..functions import *
from ..graph import *
from ..plots import *

from app import app

from .. import config

args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"AgTable",
            "name":"measure",
            "description":"Metrics of quality control to compute",
            "header":[
                { "headerName": "Name", "field":"name", "editable": True },
                { "headerName": ".var", "field":"var", "editable": True,
                  "cellEditor": "agSelectCellEditor",
                  "cellEditorParams": {"values": {"function":"list(config.adata.var.columns.values)"}},
                },
                { "headerName": "Pattern", "field":"pattern", "editable": True },
                { "headerName": "Style", "field":"style", "editable": True,
                  "cellEditor": "agSelectCellEditor",
                  "cellEditorParams": {"values": ['counts','n_expressed_genes','proportion']},
                }
            ],
            "value":[
                {
                    "name":"counts", "var":"", "pattern":"", "style":"counts",
                },
                {
                    "name":"n_expressed_genes", "var":"", "pattern":"", "style":"n_expressed_genes",
                }
            ],
            "addRows":{"name":"", "var":"", "pattern":"", "style":"counts"}
            ,
            "deleteRows": True
        },
    ],

    "postexecution" : [
        {
            "input":"Dropdown",
            "name":"batch",
            "description":"Observable to use",
            "value":None,
            "clearable":True,
            "options": {"function":"[i for i,j in zip(config.adata.obs.columns.values, config.adata.obs.dtypes) if j in ['str','object','category','int']]"}
        },
        {
            "input":"AgTable",
            "name":"thresholds",
            "description":"Thresholds to apply to the data",
            "header":[
                { "headerName": "Metric", "field":"metric", "editable": False },
                { "headerName": "Batch", "field":"batch", "editable": False },
                { "headerName": "Min", "field":"min", "editable": True },
                { "headerName": "Max", "field":"max", "editable": True },
                { "headerName": "Genes", "field":"genes", "editable": False },
            ],
            "value":{"function":"qc_data()"},
        },
        {
            "input":"Dropdown",
            "name":"plot_style",
            "description":"Chose between violin or 2D scatterplot.",
            "value":"violin",
            "clearable":False,
            "options":["violin","scatter"],
        },
        {
            "input":"Dropdown",
            "name":"plot_x",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":{"function": "qc_measures()[0]"},
            "clearable":False,
            "options":{"function": "qc_measures()"},
            "visible":{"function": "config.adata.uns[config.selected]['parameters']['plot_style'] == 'scatter'"}
        },
        {
            "input":"Dropdown",
            "name":"plot_y",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":{"function": "qc_measures()[0]"},
            "clearable":False,
            "options":{"function": "qc_measures()"},
            # "visible":{"function": "config.adata.uns[config.selected]['parameters']['plot_style'] == 'scatter'"}
        },
        {
            "input":"Dropdown",
            "name":"plot_color",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":None,
            "clearable":True,
            "options":{"function": "list(config.adata.obs.columns.values)"},
            "visible":{"function": "config.adata.uns[config.selected]['parameters']['plot_style'] == 'scatter'"}
        },
    ]

}

def qc_measures():

    selected = config.selected
    adata = config.adata

    return [f"{selected}--{i['name']}" for i in adata.uns[selected]["parameters"]["measure"]]

def qc_data():

    selected = config.selected
    adata = config.adata

    parameters = adata.uns[selected]["parameters"]

    # #Reset table
    # if "batch" in modified_arg.keys():
    #     if "batch" in config.adata.uns[config.selected]["parameters"].keys():
    #         del config.adata.uns[config.selected]["parameters"]["batch"]

    print(parameters["measure"])
    measures = [i["name"] for i in parameters["measure"]]

    data = []
    for i in parameters["measure"]:
        genes = qc_get_genes(i.copy())
        if i["name"] in measures:
            j = f"{config.selected}--{i['name']}"
            data.append({"metric":i["name"], "batch":" ","min":adata.obs[j].min(),"max":adata.obs[j].max(),"genes":genes})

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

def f_qc(adata, inputArgs, kwargs):
        
    X = inputArgs["obsm"]

    d = {"obs":{}}
    for l in kwargs["measure"]:
        genes = qc_get_genes(l)
        if genes == "all":
            x = range(X.shape[0])
        else:
            x = np.isin(adata.var[l["var"]].values, genes)

        if "counts" == l["style"]:
            d["obs"][l["name"]] = np.array(X[:,x].sum(axis=1)).reshape(-1)
        elif "n_expressed_genes" == l["style"]:
            d["obs"][l["name"]] = np.array((X[:,x] > 0).sum(axis=1)).reshape(-1)
        elif "proportion" == l["style"]:
            p = np.array(config.adata.X[:,x].sum(axis=1)).reshape(-1)/np.array(config.adata.X.sum(axis=1)).reshape(-1)
            d["obs"][l["name"]] =  p

    return d

def get_from_table(table, filter={}, properties=[]):
    properties_dict = {i:[] for i in properties}
    for i in table:
        
        add = True
        for j,k in filter.items():

            if i[j] not in k:
                add = False            

        if add:
            for j in properties_dict.keys():
                properties_dict[j].append(i[j])

    return properties_dict

def plot_qc():

    params = config.adata.uns[config.selected]["parameters"]

    if params["plot_style"] == "violin":

        if params['batch'] != None:
            x = config.adata.obs[params['batch']].values
        else:
            x = [0 for i in range(config.adata.shape[0])]

        y = config.adata.obs[params['plot_y']].values

        thresholds = get_from_table(params['thresholds'],filter={"metric":[params['plot_y'].split("--")[1]]},properties=["batch","min","max"])
        if " " in thresholds['batch']:
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
        fig.layout["yaxis"]["title"] = params['plot_y']

        return dcc.Graph(figure=fig)
        
    else:

        x = config.adata.obs[params['plot_x']].values
        y = config.adata.obs[params['plot_y']].values

        if params['plot_color'] != None:
            color = config.adata.obs[params['plot_color']].values

            fig = px.scatter(x=x,y=y,color=color)
            fig.layout["legend"]["title"]["text"] = params['plot_color']
        else:
            fig = px.scatter(x=x,y=y)

        fig.layout["xaxis"]["title"] = params['plot_x']
        fig.layout["yaxis"]["title"] = params['plot_y']

        return plot_center(dcc.Graph(figure=fig, style={"width": "60vw", "height": "50vw"}))

config.methods["qc"] = {

    "properties":{"method":"qc","type":"QC","recompute":False,"color":"blue"},

    "args": args,

    "function":f_qc,

    "plot":plot_qc,

}