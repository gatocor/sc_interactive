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

from . import config

class NpEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.integer):
            return int(obj)
        if isinstance(obj, np.floating):
            return float(obj)
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super(NpEncoder, self).default(obj)

def make_graph_path(name):
    return "./__graph__/"+name.split("/")[-1].split(".")[0]+".json"

def make_file_path(name):
    return "./__graph__/"+name.split("/")[-1].split(".")[0]+"_qc"+".h5ad"

def prevent_race(name):

    node = get_node(config.selected)
    if not node['data']['computed'] or node['data']['method'] != name:
        raise PreventUpdate()
    
def make_nodes_summaries(inplace=True):

    for node in node_names():
        if "Raw" != node:
            make_node_summary(node, inplace=True)

def make_node_summary(name, inplace=True):

    node = get_node(name)

    summary = f"{node['data']['id']}\nMethod:{node['data']['method']}\n\n"
    for prop in config.arguments[node['data']['method']]():   
        if type(prop) != str:
            if "summary" in prop.keys():
                m = node['data']['parameters'][prop['name']]
                if m == None:
                    m = '[]'
                summary += f"{prop['name']}:\n "+str(m[1:-1].replace(',','\n '))

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

def get_ancestors(name):
    l = []
    node = get_node(name)
    while node['data']['parameters']['input'] not in [None,"Raw"]:
        node = get_node(node['data']['parameters']['input'])
        l.insert(0,node)

    return l

def get_edges():
    return [i for i in config.graph if 'target' in i['data'].keys()]

def node_rm(name):

    #Remove info from node
    node = get_node(name)
    if node['data']['computed']:
        config.functions_method_rm[node['data']['method']](name)

    l = []
    for i,node in enumerate(config.graph):
        if 'id' in node['data'].keys():
            if node['data']['id'] != name:

                if name == node['data']['parameters']['input']: #Remove input from cells that have this node as input
                    node['data']['parameters']['input'] = None

                l.append(node)

        else:
            if node['data']['target'] != name and node['data']['source'] != name:
                l.append(node)

    config.graph = l

    return

def node_rename(name_old, name_new):

    for pos,i in enumerate(config.graph):
        if 'id' in i['data'].keys():
            if i['data']['parameters']['input'] == name_old: #Change all input names
                config.graph[pos]['data']['parameters']['input'] = name_new

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
        for node in get_nodes():
            if node['data']['parameters']['input'] in exclude:
                exclude.append(node['data']['id'])

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
            for j in config.arguments[i['data']['method']]():
                if j['input'] == "MeasureTable" and str2list:
                    params[j['name']] = params[j['name']][1:-1].split(",")
                    if params[j['name']][0] == '':
                        params[j['name']] = []

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
    config.graph[pos]['data']['parameters'] = config.active_node_parameters.copy()

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
            if name == node['data']['parameters']['input']:
                deactivate_node(node['data']['id'])
                node = get_node(node['data']['id'])
                if node['data']['computed']:
                    config.functions_method_rm[node['data']['method']](node['data']['id'])
                deactivate_downstream(node['data']['id'])

def image_unselected(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'.png'
    
def image_selected(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'_selected.png'

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

        return [color_map[i] for i in labels]

    elif type(labels[0]) in [str,object,int,bool,np.bool_,np.int_,"category"]:

        color_sequence = px.colors.qualitative.Plotly

        # Create a dictionary to map categories to colors
        color_map = {
            category: color_sequence[i % len(color_sequence)]
            for i, category in enumerate(np.unique(labels))
        }

        return [color_map[i] for i in labels]
    
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