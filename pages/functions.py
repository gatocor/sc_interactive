import scanpy as sc
import numpy as np
import pandas as pd
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html
import dash_daq as daq
from dash import dash_table
import dash_renderjson
import plotly.graph_objs as go
import plotly.express as px
from dash.exceptions import PreventUpdate
import json
import os

from . import config
class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        if isinstance(obj, pd.Categorical):
            return obj.tolist()
        
        return super(NpEncoder, self).default(obj)

def make_folders_structure(name):

    basedir = f"../{name}"
    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    adatadir = f"../{name}/h5ad"
    if not os.path.isdir(adatadir):
        os.mkdir(adatadir)

    reportdir = f"../{name}/report"
    if not os.path.isdir(reportdir):
        os.mkdir(reportdir)

    reportimagesdir = f"../{name}/report/figures"
    if not os.path.isdir(reportimagesdir):
        os.mkdir(reportimagesdir)

    return

def get_figures(fig):
    l = []
    try:
        f = go.Figure(fig)
        l.append(f)
    except:
        if type(fig) == list:
            for i in fig:
                l += get_figures(i)
        elif type(fig) == dict:
            for j,k in fig.items():
                l += get_figures(k)

    return l

def make_graph_path(name):
    return "./__graph__/"+name.split("/")[-1].split(".")[0]+".json"

def make_file_path(name):
    return "./__graph__/"+name.split("/")[-1].split(".")[0]+"_qc"+".h5ad"

def hist_vline(data,bins):

    hist_values, hist_bins = np.histogram(data, bins=bins)
    tallest_bin_height = np.max(hist_values)

    return tallest_bin_height

def qualitative_colors(labels):

    if type(labels[0]) not in [str,object,int,bool,np.bool_,np.int_,"category"]:

        return labels

    elif type(labels[0]) in [int,np.int_] and len(labels)>100:

        return labels
    
    elif type(labels[0]) in [int,np.int_]:

        color_map = {
            True: "lightblue", False: "red"
        }

        return np.vectorize(color_map.get)(labels)

    elif type(labels[0]) in [str,object,int,bool,np.bool_,np.int_,"category"]:

        color_sequence = px.colors.qualitative.Plotly

        # Create a dictionary to map categories to colors
        color_map = {
            category: color_sequence[i % len(color_sequence)]
            for i, category in enumerate(np.unique(labels))
        }

        return np.vectorize(color_map.get)(labels)
    
def get_batch_keys():
    return [i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in [str,object,int,"category"]]

def make_thresholds_table(threshold_measures,batches,cols_add):
    
    l = []
    for i in threshold_measures:
        for j in cols_add:
            l.append(f"{i} {j}")

    return make_table(l,batches)

def make_table(col_names,row_names):
    cols = [{"name":"Batch", "id":"Batch", "deletable": False, "editable": False}]+[{"name":i, "id":i, "deletable": False, "editable": True} for i in col_names]
    data = []
    for j in row_names:
        d = {"Batch":j}
        for i in col_names:
            d[i]=0
        data.append(d)

    return cols,data

def set_table_value(data, row_name, col_name, val):
    for i,row in enumerate(data):
        if row['Batch'] == row_name:
            data[i][col_name] = val

def set_table_column(data, col_name, val):
    for i,row in enumerate(data):
        data[i][col_name] = val

def get_table_value(data, row_name, col_name):
    for i,row in enumerate(data):
        if row['Batch'] == row_name:
            return data[i][col_name]

def get_table_column(data, col_name):
    return [i[col_name] for i in data]

def json_serializable(uns):

    d = "{"

    for i,object in uns.items():
        if type(object) == dict:
            d = f"{d}'{i}':{json_serializable(object)},"
        elif type(object) in [int, float, list, np.int_, np.float_, bool]:
            n = str(type(object)).split("class '")[-1].split("'")[0]
            d = f"{d}'{i}':{object},"
        elif type(object) in [str, list]:
            if i == "summary":
                o = object.split('\n')[0]
                d = f"{d}'{i}':'{o}',"
            elif i == "image":
                o = object
                if '_selected.png' in o:
                    o = object.split('_selected.png')[0]+'.png'
                d = f"{d}'{i}':'{o}',"
            elif i != "summary":
                n = str(type(object)).split("class '")[-1].split("'")[0]
                d = f"{d}'{i}':'{object}',"
        elif object == None:
            d = f"{d}'{i}':None,"
        else:
            n = str(type(object)).split("class '")[-1].split("'")[0]
            d = f"{d}'{i}':'{object}',"

    d = d[:-1]+"}"

    return d

def f_qc_table_pattern(adata):

    return [{"Concept":i, "Pattern":str(adata.uns["gene_lists"][i]["Pattern"]),"Genes":str(adata.uns["gene_lists"][i]["Genes"])} for i in adata.uns["gene_lists"]]

def make_adata_layout(adata):
    layout=[]
    if type(adata) != type(None):
        if "GeneNamesKey" in adata.uns.keys():
            layout = [
                        dbc.Row(
                [
                    dbc.Col(html.H1(".X"), width="auto"),
                ],
                justify="left",
                className="mb-4"
            ),
            dbc.Row(
                [
                    dash_table.DataTable(
                        id='table_x',
                        columns=[
                            {"name": i, "id": i, "deletable": False, "editable": False} for i in ["cells","genes","dtype"]
                        ],
                        data=[{"cells": str(adata.X.shape[0]), "genes": str(adata.X.shape[1]), "dtype": str(type(adata.X))}],
                        editable=False,
                        row_deletable=False,
                        fixed_rows={'headers': True},
                        style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                        style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                    ),
                ],
                style={"margin-bottom":"1cm"}
            ),
            dbc.Row(
                [
                    dbc.Col(html.H1(".obs"), width="auto"),
                ],
                justify="left",
                className="mb-4"
            ),
            dbc.Row(
                [
                    dash_table.DataTable(
                        id='table_obs',
                        columns=[
                            {"name": i, "id": i, "deletable": False, "editable": False} for i in adata.obs.columns.values
                        ],
                        data=adata.obs.to_dict("records"),
                        editable=True,
                        row_deletable=False,
                        fixed_rows={'headers': True},
                        style_table={'overflowY': 'auto', 'height': '500px', 'overflowX': 'auto'},
                        style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(html.H1(".var"), width="auto"),
                ],
                justify="left",
                className="mb-4"
            ),
            dbc.Row(
                [
                    dash_table.DataTable(
                        id='table_var',
                        columns=[
                            {"name": i, "id": i, "deletable": False, "editable": False} for i in adata.var.columns.values
                        ],
                        data=adata.var.sort_values(adata.uns["GeneNamesKey"]).to_dict("records"),
                        editable=True,
                        row_deletable=False,
                        fixed_rows={'headers': True},
                        style_table={'overflowY': 'auto', 'height': '500px', 'overflowX': 'auto'},
                        style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                    ),
                ]
            ),
            dbc.Row(
                [
                    dbc.Col(html.P("Gene lists."), width="auto"),
                ],
                justify="left",
                className="mb-4"
            ),
            dbc.Row(
            dash_table.DataTable(
                    id='qc_pattern_table',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": True} for i in ["Concept","Pattern","Genes"]
                    ],
                    data=f_qc_table_pattern(adata),
                    editable=True,
                    row_deletable=True,
                    # row_selectable="multi",
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                ),
            ),
            dbc.Row(
                dbc.Button('Add Row', id='add-row-button', n_clicks=0,
                        style={
                            "background-color":"#343A40",
                            }
                            ),
                style={"margin-bottom":"1cm"}
            ),
            dbc.Row(
                [
                    dbc.Col(html.H1(".obsm"), width="auto"),
                ],
                justify="left",
                className="mb-4"
            ),
            dbc.Row(
                [
                    dash_table.DataTable(
                        id='table_obsm',
                        columns=[
                            {"name": i, "id": i, "deletable": False, "editable": False} for i in ["name","cells","vars","dtype"]
                        ],
                        data=[{"name": i, "cells": str(adata.obsm[i].shape[0]), "vars": str(adata.obsm[i].shape[1]), "dtype": str(type(adata.obsm[i]))}
                            for i in adata.obsm
                            ],
                        editable=False,
                        row_deletable=False,
                        fixed_rows={'headers': True},
                        style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                        style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                    ),
                ],
                style={"margin-bottom":"1cm"}
            ),
            dbc.Row(
                [
                    dbc.Col(html.H1(".uns"), width="auto"),
                ],
                justify="left",
                className="mb-4"
            ),
            dbc.Row(
                dash_renderjson.DashRenderjson(id="output_uns", data=json_serializable(adata.uns), max_depth=1, invert_theme=True)
            ),
            dbc.Row(
                [
                    dbc.Button(id='gene-list-save-button', n_clicks=0, children="Save",
                                size="lg",
                                style={
                                    "background-color":"#343A40",
                                    'width': '280px', 
                                }      
                                )
                ],
                justify="center",
                className="mb-4"
            ),
            ]
    
    return layout

def de2array(de, n_genes=2):

    clusters = de["scores"].dtype.names
    labels_x = []
    for cluster in clusters:
        order = np.argsort(-np.abs(de["scores"][cluster]))
        labels_x = np.append(labels_x, de["names"][cluster][order[:n_genes]])

    labels_y = clusters

    data_array = np.zeros([n_genes*len(clusters),len(clusters)])
    color = []
    for i,cluster in enumerate(clusters):
        for j,gene in enumerate(labels_x):
            g = gene == de["names"][cluster]
            data_array[j,i] = de["scores"][cluster][g][0]
            color.append(de["scores"][cluster][g][0])

    return data_array, labels_x, labels_y

def de_markers2array(de, genes, var):

    clusters = de["scores"].dtype.names

    data_array = np.zeros([len(genes),len(clusters)])
    color = []
    for i,cluster in enumerate(clusters):
        names = de["names"][cluster]
        names = config.adata.var.loc[names,var].values
        for j,gene in enumerate(genes):
            try:
                g = gene == names
                data_array[j,i] = de["scores"][cluster][g][0]
                color.append(de["scores"][cluster][g][0])
            except:
                None

    return data_array, genes, clusters