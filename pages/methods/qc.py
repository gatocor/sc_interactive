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
            "input":"QCTable",
            "name":"measure",
            "description":"Metrics of quality control to compute",
            "value":[
                {
                    "name":"counts",
                    "var":" ",
                    "pattern":" ",
                    "style":"counts",
                    "genes":"all",
                },
                {
                    "name":"n_expressed_genes",
                    "var":" ",
                    "pattern":" ",
                    "style":"n_expressed_genes",
                    "genes":"all",
                }
            ],
        },
    ],

    "postexecution" : [
        ARGBATCH,
        {
            "input":"AgTable",
            "name":"thresholds",
            "description":"Thresholds to apply to the data",
            "header":[
                { "headerName": "Metric", "field":"metric", "editable": False },
                { "headerName": "Cluster", "field":"cluster", "editable": False },
                { "headerName": "Min", "field":"min", "editable": True },
                { "headerName": "Max", "field":"max", "editable": True },
                { "headerName": "Plot resolution", "field":"resolution", "editable": True },
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
            "value":{"function": "config.adata.obs.columns.values[0]"},
            "clearable":False,
            "options":{"function": "config.adata.obs.columns.values"},
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

    a = adata.uns[selected]["parameters"]

    # #Reset table
    # if "batch" in modified_arg.keys():
    #     if "batch" in config.adata.uns[config.selected]["parameters"].keys():
    #         del config.adata.uns[config.selected]["parameters"]["batch"]

    measures = [i["name"] for i in adata.uns[selected]["parameters"]["measure"]]

    data = []
    for i in a["measure"]:
        if i["name"] in measures:
            j = f"{config.selected}--{i['name']}"
            data.append({"metric":i["name"], "cluster":None,"min":adata.obs[j].min(),"max":adata.obs[j].max(),"resolution":adata.obs[j].std()})

    return data

def f_qc(adata, inputArgs, kwargs):
        
    X = inputArgs["obsm"]

    d = {"obs":{}}
    for l in kwargs["measure"]:
        if l["genes"] == "all":
            x = range(X.shape[0])
        else:
            x = np.isin(adata.var[l["var"]].values, l["genes"])

        if "counts" == l["style"]:
            d["obs"][l["name"]] = np.array(X[:,x].sum(axis=1)).reshape(-1)
        elif "n_expressed_genes" == l["style"]:
            d["obs"][l["name"]] = np.array((X[:,x] > 0).sum(axis=1)).reshape(-1)
        elif "proportion" == l["style"]:
            p = np.array(config.adata.X[:,x].sum(axis=1)).reshape(-1)/np.array(config.adata.X.sum(axis=1)).reshape(-1)
            d["obs"][l["name"]] =  p

    return d

def plot_qc():

    params = config.adata.uns[config.selected]["parameters"]

    if params["plot_style"] == "violin":

        if params['batch']:
            x = config.adata.obs[params['batch']].values
        else:
            x = [0 for i in range(config.adata.shape[0])]

        y = config.adata.obs[params['plot_y']].values

        l = make_Graph([plot_violin(x,y)], x_label=params['batch'], y_label=params['plot_y'])

        return l
        
    else:

        x = config.adata.obs[params['plot_x']].values
        y = config.adata.obs[params['plot_y']].values
        color = config.adata.obs[params['plot_color']].values

        l = make_Graph([plot_scatter(x, y, color)], x_label=params['plot_x'], y_label=params['plot_y'])

        return l

config.methods["qc"] = {

    "properties":{"method":"qc","type":"QC","recompute":False},

    "args": args,

    "function":f_qc,

    "plot":plot_qc,
}