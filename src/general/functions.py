import scanpy as sc
from scipy.sparse import issparse
import numpy as np
import pandas as pd
import json
import os
from shutil import rmtree, copyfile
from copy import deepcopy

from dash import dash_table
from dash import dcc
from dash import html
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_daq as daq
import dash_ag_grid as dag
import dash_renderjson

import plotly.graph_objs as go
import plotly.express as px
from plotly.subplots import make_subplots

from . import config
from .constants import *

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

def args_color():

    l = []
    for i in config.adata.obs.columns.values:
        l.append(f"obs.{i}")
    for i in config.adata.var.columns.values:
        l.append(f"var.{i}")
    for i in config.adata.obsm.keys():
        l.append(f"obsm.{i}")

    return l

def args_color_var():

    if get_node(config.selected)['data']['plot']['color'] == None:
        return [None]
    elif get_node(config.selected)['data']['plot']['color'].startswith("var."): 
        var = get_node(config.selected)['data']['plot']['color'][4:]
        return config.adata.var[var].values
    else:
        return [None]
    
def args_color_obsm():

    if get_node(config.selected)['data']['plot']['color'] == None:
        return [None]
    elif get_node(config.selected)['data']['plot']['color'].startswith("obsm."):
        obsm = get_node(config.selected)['data']['plot']['color'][5:]
        return np.arange(config.adata.obsm[obsm].shape[0])
    else:
        return [None]

def get_color():

    plot_args = get_node(config.selected)['data']['plot']

    if plot_args['color'] == None:
        return None
    
    elif plot_args['color'].startswith("var."):

        key = plot_args['color'][4:]
        pos = np.where(config.adata.var[key].values == plot_args["color_var"])[0]

        if plot_args["color_layer"] == "X":
            X = config.adata.X
        else:
            X = config.adata.layer[plot_args["color_layer"]]

        if issparse(X):
            return np.array(X[:,pos].todense()[:,0]).reshape(-1)
        else:
            return X[:,pos].reshape(-1)
        
    elif plot_args['color'].startswith("obsm."):
        key = plot_args['color'][5:]
        return config.adata.obsm[key][:,plot_args['color_obsm_dimension']]

    elif plot_args['color'].startswith("obs."):
        key = plot_args['color'][4:]
        return np.array(config.adata.obs[key].values)

def get_representation_styles():
    if config.adata.obsm[get_node(config.selected)['data']['plot']['plot_representation']].shape[1]>2: 
        return ['scatter_2d','scatter_3d']
    else: 
        return ['scatter_2d']

def get_representation(color=None):

    plot_args = get_node(config.selected)["data"]["plot"]

    X = config.adata.obsm[plot_args["plot_representation"]]

    if plot_args["plot_style"] == "scatter_2d":

        fig = make_subplots(rows=plot_args['plot_n_components']-1, cols=plot_args['plot_n_components']-1, shared_yaxes=True, shared_xaxes=True)
        for i in range(plot_args['plot_n_components']-1):

            for j in range(i+1,plot_args['plot_n_components']):
                x = X[:,i]
                y = X[:,j]
                if plot_args['plot_n_components'] == 2:
                    x = X[:,plot_args["plot_dimension_x"]]
                    y = X[:,plot_args["plot_dimension_y"]]

                fig.add_traces(
                        list(px.scatter(
                                    x=x,
                                    y=y,
                                    color=color,
                                    labels=None
                        ).select_traces()),
                        rows=j, cols=i+1
                )

    elif plot_args["plot_style"] == "scatter_3d":
        x = X[:,plot_args["plot_dimension_x"]]
        y = X[:,plot_args["plot_dimension_y"]]
        z = X[:,plot_args["plot_dimension_z"]]

        fig = px.scatter_3d(
                    x=x,
                    y=y,
                    z=z,
                    color=color,
                )

    return fig

def matrix_options():

    l = ["X"]
    if len(config.adata.layers.keys()):
        l.append("layers")
    if len(config.adata.obsm.keys()):
        l.append("obsm")

    return l

def matrix_keys():

    l = []
    if config.active_node_parameters['matrix'] == "layers":
        l = [i for i in config.adata.layers.keys()]
    elif config.active_node_parameters['matrix'] == "obsm":
        l = [i for i in config.adata.obsm.keys()]

    if l != []:
        return l
    else:
        return [None]

def get_matrix(matrix,key):
    if matrix == "X":
        l = config.adata.X
    elif matrix == "layers":
        l = config.adata.layers[key]
    elif matrix == "obsm":
        l = config.adata.layers[key]

    return l

def name_analysis_folder(file):
    
    if file.endswith(ENDH5AD):
        file = file[:-len(ENDH5AD)]

    if file.endswith(ENDANALYSIS):
        file = file[:-len(ENDANALYSIS)]

    return file + ENDANALYSIS

def create_analysis(basedir, raw):

    if not os.path.isdir(basedir):
        os.mkdir(basedir)

    adatadir = f"{basedir}/h5ad"
    if not os.path.isdir(adatadir):
        os.mkdir(adatadir)

    reportdir = f"{basedir}/report"
    if not os.path.isdir(reportdir):
        os.mkdir(reportdir)

    reportimagesdir = f"{basedir}/report/figures"
    if not os.path.isdir(reportimagesdir):
        os.mkdir(reportimagesdir)

    if raw:
        dst = f"{basedir}/h5ad/Raw.h5ad"
        copyfile(raw, dst)

    file = f"{basedir}/analysis.json"
    with open(file,"w") as outfile:
        json_object = json.dumps([RAWNODE], indent=4, cls=NpEncoder)
        outfile.write(json_object)

    file = f"{basedir}/report/report.md"
    with open(file,"w") as outfile:
        outfile.write("")

    return

def load_analysis(folder):
    
    graphfile = f"{folder}/analysis.json"
    with open(graphfile,"r") as outfile:
        config.graph = json.load(outfile)

    h5adfile = f"{folder}/h5ad/Raw.h5ad"
    load_adata(h5adfile)

    file = f"{folder}/report/report.md"
    with open(file,"r") as outfile:
        config.report = outfile.read()

    config.analysis_folder = folder
    config.selected = 'Raw'

    if not config.adata.raw:
        config.adata.raw = config.adata

    if "Raw" not in config.adata.obsm.keys():
        config.adata.obsm["Raw"] = config.adata.raw.X

def get_figures(fig):
    l = []
    try:
        f = go.Figure(data=fig)
        l.append(f)
    except:
        if type(fig) == list:
            for i in fig:
                l += get_figures(i)
        elif type(fig) == dict:
            for j,k in fig.items():
                l += get_figures(k)

    return l

def deep_substitute(object, origin, target):

    if object == origin:
        return target
    
    elif type(object) == list:
        for i,j in enumerate(object):
            object[i] = deep_substitute(j, origin, target)

    try:
        for i,j in object.items():
            object[i] = deep_substitute(j, origin, target)
    except:
        None

    return object

# def adapt_adata_saving():

#     config.adata.uns = deep_substitute(config.adata.uns, origin = None, target = "__None__")

#     for j in [i["data"] for i in get_nodes()]:
#         if j["id"] in config.adata.uns.keys():
#             for k in config.methods[j["method"]]["args"]["execution"]:
#                 if k["input"] == "AgTable":
#                     config.adata.uns[j["id"]]["parameters"][k["name"]] = pd.DataFrame(config.adata.uns[j["id"]]["parameters"][k["name"]])

#             for k in config.methods[j["method"]]["args"]["postexecution"]:
#                 if k["input"] == "AgTable":
#                     config.adata.uns[j["id"]]["parameters"][k["name"]] = pd.DataFrame(config.adata.uns[j["id"]]["parameters"][k["name"]])

# def adapt_adata_loaded():

#     for j in [i["data"] for i in get_nodes()]:
#         if j["id"] in config.adata.uns.keys():
#             for k in config.methods[j["method"]]["args"]["execution"]:
#                 if k["input"] == "AgTable":
#                     config.adata.uns[j["id"]]["parameters"][k["name"]] = config.adata.uns[j["id"]]["parameters"][k["name"]].to_dict("records")

#             for k in config.methods[j["method"]]["args"]["postexecution"]:
#                 if k["input"] == "AgTable":
#                     config.adata.uns[j["id"]]["parameters"][k["name"]] = config.adata.uns[j["id"]]["parameters"][k["name"]].to_dict("records")

#     config.adata.uns = deep_substitute(config.adata.uns, origin = "__None__", target = None)

def get_name(name):
    return config.selected+"--"+name

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

def get_from_table(table, filter={}, properties={}):
    properties_dict = {i:[] for i in properties.keys()}
    for i in table:
        
        add = True
        for j,k in filter.items():

            if i[j] not in k:
                add = False            

        if add:
            for j in properties_dict.keys():
                properties_dict[j].append(properties[j](i[j]))

    return properties_dict

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

def clean_arguments(name, args, d):
    for i in args:
        if "recomputeAfter" in i.keys():
            if name in i["recomputeAfter"]:
                try:
                    del d[i['name']]
                except:
                    None

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
                dash_renderjson.DashRenderjson(id="output_uns", data=json_serializable(deepcopy(adata.uns)), max_depth=1, invert_theme=True)
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

def get_value(arg):

    if arg['input'] == 'Input':
        input = arg["properties"]["value"]
    elif arg['input'] == 'Dropdown':
        input = arg["properties"]["value"]
    elif arg['input'] == 'BooleanSwitch':
        input = arg["properties"]["on"]
    elif arg['input'] == 'AgTable':
        input = arg["properties"]["data"]

    return input

############################################################################################################################################
############################################################################################################################################
# Functions
############################################################################################################################################
############################################################################################################################################
def save_graph():

    unselect_node(config.selected)
    file = f"{config.analysis_folder}/analysis.json"
    with open(file,"w") as outfile:
        json_object = json.dumps(config.graph, indent=4, cls=NpEncoder)
        outfile.write(json_object)
    select_node(config.selected)

def save_adata():

    # adapt_adata_saving()

    pos = get_node_pos(config.selected)
    file = f"{config.analysis_folder}/h5ad/{config.graph[pos]['data']['h5ad_file']}"
    config.adata.write(file)

    # adapt_adata_loaded()

def load_adata(file):
    config.adata = sc.read(file)
    # adapt_adata_loaded()

def save_report():
    
    file = f"{config.analysis_folder}/report/report.md"
    with open(file,"w") as outfile:
        outfile.write(config.report)

def make_nodes_summaries(inplace=True):

    for node in node_names():
        if "Raw" != node:
            make_node_summary(node, inplace=True)

def make_node_summary(name, inplace=True):

    node = get_node(name)

    summary = f"{node['data']['id']}\nMethod:{node['data']['method']}\n\n"
    for prop in config.methods[node['data']['method']]["args"]():   
        if type(prop) != str:
            if "summary" in prop.keys():
                m = str(node['data']['parameters'][prop['name']])
                summary += f"{prop['name']}:\n "+str(m.replace(',','\n '))+"\n"

    if inplace:
        pos = get_node_pos(name)
        config.graph[pos]['data']['summary'] = summary
        return
    else:
        return summary

def get_node_pos(name):
    for j,i in enumerate(config.graph):
        if 'id' in i['data'].keys():
            if i['data']['id'] == name:
                return j
        
    return None

def get_nodes():
    return [i for i in config.graph if 'id' in i['data'].keys()]

def get_node_ancestors(name):
    l = []
    node = get_node(name)
    while node['data']['parameters']['input'] not in [None]:
        node = get_node(node['data']['parameters']['input'])
        l.insert(0,node)

    return l

def get_node_children(name):

    l = []
    for node in get_nodes():
        if node["data"]["parameters"]["input"] == name:
            l.append(node)

    return l

def get_edges():
    return [i for i in config.graph if 'target' in i['data'].keys()]

def node_rename(name_old, name_new):

    for pos,i in enumerate(config.graph):
        if 'id' in i['data'].keys():
            if i['data']['id'] == name_old:
                config.graph[pos]['data']['id'] = name_new
                config.graph[pos]['data']['summary'] = make_node_summary(name_new, inplace=False)
        else:
            if i['data']['target'] == name_old:
                config.graph[pos]['data']['target'] = name_new

            if i['data']['source'] == name_old:
                config.graph[pos]['data']['source'] = name_new

def node_names(exclude_downstream_from_node=None):

    exclude = []
    if exclude_downstream_from_node != None:
        exclude.append(exclude_downstream_from_node)
        for edge in get_edges():
            None
            # if edge['data']['parameters']['input'] in exclude:
            #     exclude.append(node['data']['id'])

    d = [i['data']['id'] for i in config.graph if 'id' in i['data'].keys()]
    d = [i for i in d if i not in exclude] #exclude

    return d

def node_update_pos(graph, id, pos_new):
    
    for pos,i in enumerate(graph):
        if 'id' in i['data'].keys():
            if i['data']['id'] == id:
                config.graph[pos]['position'] = pos_new

def get_node(name):
    for i in get_nodes():
        if i['data']['id'] == name:
            return i
        
    return None

def get_node_parameters(name, str2list=False):

    for i in get_nodes():
        if i['data']['id'] == name:
            params = i['data']['parameters'].copy()

            return params
        
    return None

def edge_add(source,target):

    if source != None and target != None:
        for i,val in enumerate(config.graph):
            #Change edge input
            if 'target' in val['data'].keys(): #Check is an edge
                if val['data']['source'] == source and val['data']['target'] == target : #there already an edge
                    config.graph.pop(i)
        
        config.graph.append(
            {
                'data':{'source':source,'target':target},
            }
            )

def edge_rm(source,target):

    if source != None and target != None:
        for i,val in enumerate(config.graph):
            #Change edge input
            if 'target' in val['data'].keys(): #Check is an edge
                if val['data']['source'] == source and val['data']['target'] == target : #there already an edge
                    config.graph.pop(i)

def update_node(name):
    pos = get_node_pos(name)
    config.graph[pos]['data']['parameters'] = deepcopy(config.active_node_parameters)

def deactivate_node(name):
    pos = get_node_pos(name)
    config.graph[pos]['data']['opacity'] = .3
    config.graph[pos]['data']['computed'] = False

def activate_node(name):
    pos = get_node_pos(name)
    config.graph[pos]['data']['opacity'] = 1
    config.graph[pos]['data']['computed'] = True

def deactivate_downstream(name):

    for node in get_nodes():
        if "input" in node['data']['parameters'].keys():
            if name == node['data']['parameters']['input'] and node['data']['computed']:
                deactivate_node(node['data']['id'])
                node = get_node(node['data']['id'])
                deactivate_downstream(node['data']['id'])

def unselect_node(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'.png'
    
def select_node(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'_selected.png'

def list_observables():
    l = list(config.adata.obs.columns.values)
    ancestors = get_nodes()
    for i in ancestors:
        if "obs" in i["data"].keys():
            for j in i["data"]["obs"].keys():
                l += [f"{i['data']['name']}--{j}"]

    return l

def prevent_race(name,computed=True,method=True):

    node = get_node(config.selected)
    if not node['data']['computed'] and computed:
        raise PreventUpdate()

    if node['data']['method'] != name and method:
        raise PreventUpdate()
    
def is_computed():

    return get_node(config.selected)["data"]["computed"]

def graph2table():
    return [{"Name":i['data']['id'],"Type":i['data']['type'],"Method":i['data']['method']} for i in get_nodes() if i['data']['id'] != 'Raw']

# def set_active_node_parameters(args):
#     config.active_node_parameters = args

def set_parameters(args, arg_type):

    pos = get_node_pos(config.selected)
    for i,j in args.items():
        config.graph[pos]["data"][arg_type][i] = j

    # for i,j in args.items():
    #     config.adata.uns[config.selected][arg_type][i] = j

    # config.adata.uns[config.selected]["scinteractive"] = True
    # config.adata.uns[config.selected]["method"] = get_node(config.selected)["data"]["method"]

def get_args(name):

    d = {}
    #obs
    obs_key = f"{name}--"
    l = [i for i in config.adata.obs.columns.values if i.startswith(obs_key)]   
    if len(l) > 0:   
        d["obs"] = config.adata.obs[l]

    #var
    var_key = f"{name}--"
    l = [i for i in config.adata.var.columns.values if i.startswith(var_key)]            
    if len(l) > 0:   
        d["var"] = config.adata.var[l]

    #obsm
    obsm_key = name
    if obsm_key in config.adata.obsm.keys():
        d["obsm"] = config.adata.obsm[obsm_key]
            
    #uns
    uns_key = name
    if uns_key in config.adata.uns.keys():
        d["uns"] = config.adata.uns[uns_key]

    return d

def del_adata_node(name):

    #obs
    obs_key = f"{name}--"
    l = [i for i in config.adata.obs.columns.values if i.startswith(obs_key)]   
    config.adata.obs.drop(l, axis=1, inplace=True)

    #var
    var_key = f"{name}--"
    l = [i for i in config.adata.var.columns.values if i.startswith(var_key)]            
    config.adata.var.drop(l, axis=1, inplace=True)

    #obsm
    obsm_key = name
    if obsm_key in config.adata.obsm.keys():
        del config.adata.obsm[obsm_key]
            
    #uns
    uns_key = name
    if uns_key in config.adata.uns.keys():
        del config.adata.uns[uns_key]

def is_node(node):
    return 'id' in node['data'].keys()

def is_edge(node):
    return 'source' in node['data'].keys()

def node_reassign_input(name,new_input):

    for i,node in enumerate(config.graph):
        if is_node(node):
            if name == node['data']['id']: #Remove input from cells that have this node as input
                edge_rm(config.graph[i]["data"]['parameters']['input'],name)
                edge_add(new_input,name)
                config.graph[i]["data"]['parameters']['input'] = new_input

def nodes_reassign_input(old_input,new_input):

    for i,node in enumerate(config.graph):
        if is_node(node):
            if old_input == node['data']['parameters']['input']: #Remove input from cells that have this node as input
                config.graph[i]["data"]['parameters']['input'] = new_input
                edge_rm(old_input,node["data"]["id"])
                edge_add(new_input,node["data"]["id"])

def node_rm(name):

    input = get_node(name)["data"]["parameters"]["input"]

    nodes_reassign_input(name, input)

    l = []
    for i,node in enumerate(config.graph):
        if 'id' in node['data'].keys():
            if node['data']['id'] != name:

                l.append(node)
        else:
            if node['data']['target'] != name and node['data']['source'] != name:
                l.append(node)

    # #Remove info from adata
    # del_adata_node(name)

    config.graph = l

    return

# def set_output(args):

#     #obs
#     obs_key = f"{config.selected}--"
#     if "obs" in args.keys():
#         for i,j in args["obs"].items():
#             i = f"{obs_key}{i}"
#             config.adata.obs[i] = j

#     #var
#     var_key = f"{config.selected}--"
#     if "var" in args.keys():
#         for i,j in args["var"].items():
#             i = f"{var_key}{i}"
#             config.adata.var[i] = j

#     #obsm
#     if "obsm" in args.keys():
#         config.adata.obsm[config.selected] = args["obsm"]
            
#     #uns
#     if "uns" in args.keys():
#         for i,j in args["uns"].items():
#             config.adata.uns[config.selected][i] = j

def ag_table(id, columnDefs, rowData, deleteRows=False, suppressRowTransform=False, suppressRowClickSelection=False):

    if deleteRows:
        columnDefs = [
        {
                "headerName": "",
                "cellRenderer": "DeleteButton",
                "lockPosition":'left',
                "maxWidth":35,
                "filter": False,
                'cellStyle': {'paddingRight': 0, 'paddingLeft': 0},
                "editable":False
            },
        ] + columnDefs

    fig = dag.AgGrid(
        id=id,
        rowData=rowData,
        columnDefs=columnDefs,
        defaultColDef={"editable": True},
        deleteSelectedRows=deleteRows,
        dashGridOptions={"suppressRowTransform":suppressRowTransform, "suppressRowClickSelection":suppressRowClickSelection},
        columnSize="sizeToFit",
    )

    return fig

def new_h5ad_file():

    count = 1
    file_name = f"Analysis_{count}.h5ad"
    while os.path.exists(f"{config.analysis_folder}/h5ad/{file_name}"):
        count += 1
        file_name = f"Analysis_{count}.h5ad"

    return file_name

def clean_h5ad():

    l = os.listdir(f"{config.analysis_folder}/h5ad")
    f = [i["data"]["h5ad_file"] for i in get_nodes()]
    for file in l:
        if file not in f:
            os.remove(f"{config.analysis_folder}/h5ad/{file}")
        