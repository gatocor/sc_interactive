import scanpy as sc
import numpy as np
import pandas as pd
import dash_bootstrap_components as dbc
from dash import dcc
from dash import html
import dash_daq as daq
from dash import dash_table
import dash_renderjson
from .methods.doublet_args import doublet_args
import plotly.graph_objs as go

from . import config

def make_node_summary(name, inplace=True):

    node = get_node(name)

    summary = f"{node['data']['id']}\nMethod:{node['data']['method']}\n\n"
    for prop in config.arguments[node['data']['type']]():   
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

def get_edges():
    return [i for i in config.graph if 'target' in i['data'].keys()]

def node_rm(name):

    #Remove info from node
    node = get_node(name)
    if node['data']['computed']:
        config.functions_method_rm[node['data']['type']](name)

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
            for j in config.arguments[i['data']['type']]():
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
                    config.functions_method_rm[node['data']['type']](node['data']['id'])
                deactivate_downstream(node['data']['id'])

def image_unselected(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'.png'
    
def image_selected(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'_selected.png'

# def f_qc_base(adata):
#     if "raw" not in adata.layers.keys():
#         adata.raw = adata.copy()
#     if "total_counts" not in adata.obs.columns.values:
#         adata.obs["total_counts"] = np.array(adata.raw.X.sum(axis=1)).reshape(-1)
#     if "n_genes_by_counts" not in adata.obs.columns.values:
#         adata.obs["n_genes_by_counts"] = np.array((adata.raw.X > 0).sum(axis=1)).reshape(-1)
#     if "qc" not in adata.uns.keys():
#         adata.uns["qc"] = {
#             "total_counts":{"Minimum threshold":adata.obs["total_counts"].min(), "Maximum threshold":adata.obs["total_counts"].max()},
#             "n_genes_by_counts":{"Minimum threshold":adata.obs["n_genes_by_counts"].min(), "Maximum threshold":adata.obs["n_genes_by_counts"].max()},
#         }
#     if "total_counts" not in adata.uns["qc"].keys():
#         adata.uns["qc"]["total_counts"] = {"Minimum threshold":adata.obs["total_counts"].min(), "Maximum threshold":adata.obs["total_counts"].max()}
#     if "n_genes_by_counts" not in adata.uns["qc"].keys():
#         adata.uns["qc"]["n_genes_by_counts"] = {"Minimum threshold":adata.obs["n_genes_by_counts"].min(), "Maximum threshold":adata.obs["n_genes_by_counts"].max()}
#     if "gene_lists" not in adata.uns.keys():
#         adata.uns["gene_lists"] = {}

# def f_qc(adata,metrics,data):

#     m = [i["Name"] for i in metrics]
#     l = list(adata.uns["qc"].keys())

#     if "total_counts" not in adata.obs.columns.values:
#         adata.obs["total_counts"] = np.array(adata.raw.X.sum(axis=1)).reshape(-1)
#     if "n_genes_by_counts" not in adata.obs.columns.values:
#         adata.obs["n_genes_by_counts"] = np.array((adata.raw.X > 0).sum(axis=1)).reshape(-1)
#     if "total_counts" not in m:
#         adata.uns["qc"]["total_counts"] = {"Minimum threshold":adata.obs["total_counts"].min(), "Maximum threshold":adata.obs["total_counts"].max()}
#     if "n_genes_by_counts" not in m:
#         adata.uns["qc"]["n_genes_by_counts"] = {"Minimum threshold":adata.obs["n_genes_by_counts"].min(), "Maximum threshold":adata.obs["n_genes_by_counts"].max()}
    
#     for i in l:
#         if (i not in ["total_counts","n_genes_by_counts"]) and (i not in m):
#             del adata.uns["qc"][i]

#     for i in m:
#         if (i not in ["total_counts","n_genes_by_counts"]) and (i not in l):
#             adata.obs[i] = np.array(adata.raw.X[:,adata.var[i].values].sum(axis=1)).reshape(-1)/adata.obs["total_counts"].values
#             adata.obs[i] = np.nan_to_num(adata.obs[i].values, nan=0)
#             adata.uns["qc"][i] = {"Minimum threshold":adata.obs[i].min(), "Maximum threshold":adata.obs[i].max()}
#         else:
#             for j in data:
#                 if j["Control measure"] == i:
#                     adata.uns["qc"][i]["Minimum threshold"] = float(j["Minimum threshold"])
#                     adata.uns["qc"][i]["Maximum threshold"] = float(j["Maximum threshold"])

#     return

# def plot_qc_global_hist(adata):
#     l = []
    
#     for var_selected_data in [i for i in adata.uns["qc"]]:
#         #Plot_type
#         # Create a vertical line at the specified input value
#         hist_values, hist_bins = np.histogram(adata.obs[var_selected_data].values, bins=30)
#         tallest_bin_height = np.max(hist_values)
#         var1_vertical_line_min = go.Scatter(
#             x=[adata.uns["qc"][var_selected_data]["Minimum threshold"], adata.uns["qc"][var_selected_data]["Minimum threshold"]],
#             y=[0, tallest_bin_height],
#             mode='lines',
#             name='Min threshold',
#             line=dict(color='red', dash='dash')
#         )
#         var1_vertical_line_max = go.Scatter(
#             x=[adata.uns["qc"][var_selected_data]["Maximum threshold"], adata.uns["qc"][var_selected_data]["Maximum threshold"]],
#             y=[0, tallest_bin_height],
#             mode='lines',
#             name='Max threshold',
#             line=dict(color='green', dash='dash')
#         )

#         l += [
#                 dbc.Row(
#                     [dcc.Graph(id="Holi",
#                             figure={
#                                     "data":[
#                                         go.Histogram(
#                                             x=adata.obs[var_selected_data].values,
#                                             nbinsx=30,
#                                             name='Histogram',
#                                             marker=dict(color='blue'),
#                                             opacity=0.7
#                                         ),
#                                         var1_vertical_line_min,
#                                         var1_vertical_line_max
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
        
#     return l

# def f_options(adata,motive):
#     return [i for i in adata.obs.columns.values if i.startswith(motive)]

# def f_update_patterns(adata, pattern):

#     p = adata.uns["gene_lists"]
#     for i in pattern:
#         if i["Concept"] == '' and i["Pattern"] == '' and i["Genes"] == '':
#             adata.uns["gene_lists"][""] = {"Pattern":"", "Genes":[]}
#         else:
#             if i["Concept"] != '' and i["Concept"] not in adata.uns["gene_lists"].keys():
#                 adata.uns["gene_lists"][i["Concept"]] = {"Pattern":"", "Genes":""}
#             elif i["Concept"] != '':
#                 if i["Pattern"] != '' and i["Pattern"] != None:
#                     pp = [j.startswith(i["Pattern"]) for j in adata.var[adata.uns["GeneNamesKey"]].values]
#                     adata.uns["gene_lists"][i["Concept"]] = {"Pattern":i["Pattern"], "Genes":list(adata.var[adata.uns["GeneNamesKey"]].values[pp])}
#                 elif "[" not in i["Genes"]:
#                     adata.uns["gene_lists"][i["Concept"]] = {"Pattern":i["Pattern"], "Genes":list([j for j in i["Genes"].split(" ")])}
#                 else:
#                     adata.uns["gene_lists"][i["Concept"]] = {"Pattern":i["Pattern"], "Genes":adata.uns["gene_lists"][i["Concept"]]["Genes"]}
    
#     for i in [j for j in p.keys()]:
#         if i not in [j["Concept"] for j in pattern]:
#             try:
#                 del p[i]
#                 adata.var.drop(i,axis=1,inplace=True)
#             except:
#                 None
#         elif i != '':
#             adata.var[i] = [j in adata.uns["gene_lists"][i]["Genes"] for j in adata.var[adata.uns["GeneNamesKey"]].values]

#     return

# def f_qc_table_metrics(adata):
    
#     l = [{"Name":i} for i in adata.uns["qc"]]

#     return l

# def f_qc_table_threshold(adata):
    
#     l = [{"Control measure":i,
#           "Minimum threshold":str(adata.uns["qc"][i]["Minimum threshold"]),
#           "Maximum threshold":str(adata.uns["qc"][i]["Maximum threshold"]),
#           } for i in adata.uns["qc"]]

#     return l

# def f_qc_summary(adata):
#     statistics = {}

#     return pd.DataFrame(statistics)

def json_serializable(uns):

    d = {}

    for i in uns:
        if type(uns[i]) == dict:
            d[i] = json_serializable(uns[i])
        elif type(uns[i]) in [int, float, str, list, np.int_, np.float_]:
            d[i] = uns[i]
        else:
            d[i] = str(type(uns[i]))

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

# def make_qc_per_condition(adata):

#     if "per_condition" not in adata.uns["qc"]["total_counts"].keys():

#         return [dbc.Button(id='qc_per_condition-button', n_clicks=0, children="Add per conditon analysis",
#                         size="lg",
#                         style={
#                             "background-color":"gray",
#                             # 'width': '680px', 
#                         }      
#                         )]
#     else:
        
#         return  [
#                 dbc.Row(
#                     [
#                         dbc.Col(
#                             html.H1("Per condition Analysis")
#                         ),
#                         dbc.Col(
#                             dbc.Button(id='qc_per_condition-button', n_clicks=1, children="Remove per condition analysis",
#                                     size="lg",
#                                     style={
#                                         "background-color":"#343A40",
#                                         'width': '280px', 
#                                     }      
#                                     ),
#                             width=2
#                         ),
#                     ],
#                     style={"margin-bottom":"1cm"}
#                 ),
#                 dbc.Row(
#                         dash_table.DataTable(
#                             id='table_qc_per_condition_metrics',
#                             columns=[
#                                 {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Name"]
#                             ],
#                             data=[{"Name":i} for i in adata.uns["qc"]["total_counts"]["per_condition"].keys()],
#                             editable=False,
#                             row_deletable=True,
#                             style_table={'overflowY': 'auto', 'overflowX': 'auto'},
#                             style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
#                         ),
#                         # style={"margin-bottom":"1cm"}        
#                     ),
#                 dbc.Row(
#                     [
#                         dbc.Col(dcc.Dropdown(
#                             id = "dropdown_add_per_condition_metrics",
#                             value = None,
#                             options = [str(i) for i in adata.obs.columns.values if (adata.obs.dtypes[i] in ["category" ,object, str, int])]
#                         )),
#                         dbc.Col(dbc.Button("Add condition", id="add_qc_per_condition-button")),
#                     ],
#                     justify="Left",
#                     className="mb-4",
#                 ),
#                 dbc.Row(
#                     id="per_condition_plot",
#                     children=qc_per_condition_plots(config.adata),
#                 ),
#                 dbc.Row(
#                     dash_table.DataTable(
#                             id='table_qc_per_condition',
#                             columns=[
#                                 {"name": str(i), "id": str(i), "deletable": False, "editable": j} for i,j in zip(["Condition","Variable","Condition type","Min","Max"],[False,False,False,True,True])
#                             ],
#                             data=qc_per_condition_table(config.adata),
#                             editable=True,
#                             row_deletable=False,
#                             style_table={'overflowY': 'auto', 'overflowX': 'auto'},
#                             style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
#                         ),                    
#                 )
#             ]
    
# def qc_per_condition_plots(adata):

#     lp = []
#     for condition in adata.uns["qc"]["total_counts"]["per_condition"].keys():

#         lp += [
#                 dbc.Row(
#                     html.H1(condition),
#                     justify="left"
#                 )
#         ]

#         for var_selected_data in [i for i in adata.uns["qc"]]:                
#             #Plot_type
#             # Create a vertical line at the specified input value
#             var1_vertical_line_min = go.Scatter(
#                 x=adata.obs[condition].values,
#                 y=[adata.uns["qc"][var_selected_data]["Minimum threshold"] for i in range(len(adata.obs[var_selected_data].values))],
#                 mode='lines',
#                 name='Global min threshold',
#                 line=dict(color='red')
#             )
#             var1_vertical_line_max = go.Scatter(
#                 x=adata.obs[condition].values,
#                 y=[adata.uns["qc"][var_selected_data]["Maximum threshold"] for i in range(len(adata.obs[var_selected_data].values))],
#                 mode='lines',
#                 name='Global max threshold',
#                 line=dict(color='green')
#             )

#             custom_marker = {
#                 'symbol': 'line-ns',  # Symbol code for a horizontal line (short dash)
#                 'size': 3,  # Length of the horizontal line
#                 'color': 'blue',  # Color of the line
#                 'line': {'width': 20}  # Line width
#             }

#             var1_vertical_line_per_condition_min = go.Scatter(
#                 x=np.unique(adata.obs[condition].values),
#                 y=[adata.uns["qc"][var_selected_data]["per_condition"][condition][i]["Min"] for i in np.unique(adata.obs[condition].values)],
#                 mode='markers',
#                 marker=custom_marker,
#                 name='Local min threshold',
#                 line=dict(color='darkred')
#             )

#             var1_vertical_line_per_condition_max = go.Scatter(
#                 x=np.unique(adata.obs[condition].values),
#                 y=[adata.uns["qc"][var_selected_data]["per_condition"][condition][i]["Max"] for i in np.unique(adata.obs[condition].values)],
#                 mode='markers',
#                 name='Local max threshold',
#                 line=dict(color='darkgreen')
#             )

#             lp += [
#                     dbc.Row(
#                         [
#                             dcc.Graph(id="Holi",
#                                 figure={
#                                         "data":[
#                                             go.Violin(
#                                                 x=adata.obs[condition].values,
#                                                 y=adata.obs[var_selected_data].values,
#                                                 name='Violin',
#                                                 # marker=dict(color='blue'),
#                                                 opacity=0.7
#                                             ),
#                                             var1_vertical_line_min,
#                                             var1_vertical_line_max,
#                                             var1_vertical_line_per_condition_min,
#                                             var1_vertical_line_per_condition_max
#                                         ],
#                                         "layout":{
#                                                 'title': f'{var_selected_data}',
#                                                 'xaxis': {'title': condition},
#                                                 'yaxis': {'title': 'Count'},
#                                                 'barmode': 'overlay',
#                                                 'width':1500,
#                                                 'height':400,
#                                         }
#                                     }
#                             )
#                         ],
#                         # justify="center",
#                         # style={'width': '90%', 'margin': 'auto'}
#                     )
#                 ]
            
#     return lp

# def qc_per_condition_table(adata):
#     datas = []
#     for condition in adata.uns["qc"]["total_counts"]["per_condition"].keys():

#         for var_selected_data in [i for i in config.adata.uns["qc"]]:                
            
#             for i in np.unique(config.adata.obs[condition].values):

#                 datas.append({"Condition":condition,
#                              "Variable":var_selected_data,
#                              "Condition type":i,
#                              "Min":config.adata.uns["qc"][var_selected_data]["per_condition"][condition][i]["Min"],
#                              "Max":config.adata.uns["qc"][var_selected_data]["per_condition"][condition][i]["Max"]})

#     return datas

# def make_qc_doublet_plots(adata):

#     if adata.uns["doublets"] != {}:

#         hist_values, _ = np.histogram(adata.obs["doublet_score"], bins=30)
#         tallest_bin_height = np.max(hist_values)
#         hist_values, _ = np.histogram(adata.uns["doublets"]["doublet_score_sim"], bins=30)
#         tallest_bin_height_2 = np.max(hist_values)
#         tallest_bin_height = max(tallest_bin_height,tallest_bin_height_2)

#         plot_hist = [
#             go.Histogram(
#                 x=adata.obs["doublet_score"],
#                 nbinsx=30,
#                 name='Obs',
#                 marker=dict(color='blue'),
#                 opacity=0.7
#             ),
#             go.Histogram(
#                 x=adata.uns["doublets"]["doublet_score_sim"],
#                 nbinsx=30,
#                 name='Simulated',
#                 marker=dict(color='green'),
#                 opacity=0.7
#             ),
#             go.Scatter(
#                     x=[adata.uns["doublets"]["threshold"]["full"], adata.uns["doublets"]["threshold"]["full"]],
#                     y=[0,tallest_bin_height],
#                     mode='lines',
#                     marker={"color":"red"},
#                     name='Threshold',
#             )
#         ]

#         plot_umap = [
#                 go.Scatter(
#                     x=adata.uns["doublets"]["X_umap"][:,0],
#                     y=adata.uns["doublets"]["X_umap"][:,1],
#                     mode='markers',
#                     marker={"color":adata.obs["doublet_score"]},
#                     name='Min threshold',
#             ),
#         ]

#         return dbc.Row([
#             dbc.Col(
#                 dcc.Graph(figure={"data":plot_hist})
#             ),
#             dbc.Col(
#                 dcc.Graph(figure={"data":plot_umap})
#             )
#             ]
#         )
#     else:
#         return []

# def make_qc_doublets(adata):

#     if "doublets" not in adata.uns.keys():
#         return [dbc.Button(id='doublets-button', n_clicks=0, children="Add Doublet Analysis",
#                             size="lg",
#                             style={
#                                 "background-color":"gray",
#                             }      
#                             )]
    
#     else:
#         l = make_arguments("doublets",doublet_args(adata),{})

#         return  [
#                     dbc.Row(
#                         [
#                             dbc.Col(
#                                 html.H1("Doublet Analysis")
#                             ),
#                             dbc.Col(
#                                 dbc.Button(id='doublets-button', n_clicks=1, children="Remove Doublet Analysis",
#                                         size="lg",
#                                         style={
#                                             "background-color":"#343A40",
#                                             'width': '280px', 
#                                         }      
#                                         ),
#                                 width=2
#                             ),
#                         ],
#                         style={"margin-bottom":"1cm"}
#                     ),
#                     dbc.Col(l,width=4,
#                             style={"background-color":"lightgray"}
#                             ),
#                     dbc.Col([
#                         dbc.Row(id="qc_plot_doublets",children=make_qc_doublet_plots(adata))
#                     ]
#                     ),
#                     dash_table.DataTable(
#                             id='table_doublets_threshold',
#                             columns=[
#                                 {"name": str(i), "id": str(i), "deletable": False, "editable": j} for i,j in zip(["Condition","Threshold"],[False,True])
#                             ],
#                             data=[],
#                             editable=True,
#                             row_deletable=False,
#                             style_table={'overflowY': 'auto', 'overflowX': 'auto'},
#                             style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
#                         ),                    
#                 ]

# def qc_summary(adata):

#     summary = []

#     #qc
#     if "per_condition" not in adata.uns["qc"]["total_counts"].keys(): #global
#         for i,j in  adata.uns["qc"].items():
#             rem = adata.obs[i].values < j["Minimum threshold"]
#             rem += adata.obs[i].values > j["Maximum threshold"]
#             rem >= 0

#             summary.append({"Concept":i,"Removed cells":np.sum(rem),"Removed cells (%)":np.round(np.sum(rem)*100/len(rem),2)})
#     else: #local
#         for i,qcs in  adata.uns["qc"].items():
#             rem = np.zeros(adata.shape[0])
#             for j,cond in qcs["per_condition"].items():
#                 for k,key in  cond.items():
#                     r = (adata.obs[i].values < key["Min"]) * (adata.obs[j] == k)
#                     rem[r > 0] += 1
#                     r = (adata.obs[i].values > key["Max"]) * (adata.obs[j] == k)
#                     rem[r > 0] += 1

#                 rem >= 0

#             summary.append({"Concept":i,"Removed cells":np.sum(rem),"Removed cells (%)":np.round(np.sum(rem)*100/len(rem),2)})

#     #doublets
#     if "doublets" in adata.uns.keys() and adata.uns["doublets"] != {}:

#         if "full" in adata.uns["doublets"]["batch_key"]:
#             rem = adata.obs["doublet_score"].values > adata.uns["doublets"]["threshold"]["full"]

#             summary.append({"Concept":"doublets","Removed cells":np.sum(rem),"Removed cells (%)":np.round(np.sum(rem)*100/len(rem),2)})
#         else:
#             rem = np.zeros(adata.shape[0])
#             for i,j in adata.uns["doublets"]["threshold"].items():
#                 rem = adata.obs["doublet_score"].values > adata.uns["doublets"]["threshold"]["full"]

#             summary.append({"Concept":"doublets","Removed cells":np.sum(rem),"Removed cells (%)":np.round(np.sum(rem)*100/len(rem),2)})


#     return summary

# def qc_limit(adata):

#     rem_total = np.zeros(adata.shape[0])

#     #qc
#     if "per_condition" not in adata.uns["qc"]["total_counts"].keys(): #global
#         for i,j in  adata.uns["qc"].items():
#             rem = adata.obs[i].values < j["Minimum threshold"]
#             rem += adata.obs[i].values > j["Maximum threshold"]
#             rem >= 0

#             rem_total += rem 
#     else: #local
#         for i,qcs in  adata.uns["qc"].items():
#             rem = np.zeros(adata.shape[0])
#             for j,cond in qcs["per_condition"].items():
#                 for k,key in  cond.items():
#                     r = (adata.obs[i].values < key["Min"]) * (adata.obs[j] == k)
#                     rem[r > 0] += 1
#                     r = (adata.obs[i].values > key["Max"]) * (adata.obs[j] == k)
#                     rem[r > 0] += 1

#                 rem >= 0

#             rem_total += rem 

#     #doublets
#     if "doublets" in adata.uns.keys() and adata.uns["doublets"] != {}:
#         if "full" in adata.uns["doublets"]["batch_key"]:
#             rem = adata.obs["doublet_score"].values > adata.uns["doublets"]["threshold"]["full"]

#             rem_total += rem 
#         else:
#             rem = np.zeros(adata.shape[0])
#             for i,j in adata.uns["doublets"]["threshold"].items():
#                 rem = adata.obs["doublet_score"].values > adata.uns["doublets"]["threshold"]["full"]

#                 rem_total += rem 

#     return rem_total == 0

# def first_word(l):
#     ll = l.split(" ")
#     for i in ll:
#         if i != "":
#             return i
        
#     return l

# def get_parameter_descriptions(func):

#     docstring = func.__doc__
#     if not docstring:
#         return {}
    
#     parameters = inspect.signature(func).parameters
#     param_names = list(parameters.keys())  

#     param_descriptions = {}
#     lines = docstring.split('\n')

#     name = None
#     desc = ""
#     for line in lines:
#         l = first_word(line)
#         if any([l == i for i in param_names]):
#             if name != None:
#                 param_descriptions[name]=desc
#             l = first_word(line)
#             name = l
#             desc = ""
#         elif line.startswith("        "):
#             desc += line[8:] + "\n"
#         else:
#             if name != None:
#                 param_descriptions[name]=desc
#             name = None
    
#     return param_descriptions

# def extract_keyword_arguments_with_defaults(func):
#     signature = inspect.signature(func)
#     keyword_arguments_with_defaults = {}

#     for param_name, param in signature.parameters.items():
#         if param.default != inspect.Parameter.empty:
#             keyword_arguments_with_defaults[param_name] = param.default

#     return keyword_arguments_with_defaults
