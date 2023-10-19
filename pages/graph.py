import os
import dash
import dash_bootstrap_components as dbc
from dash import dcc, ctx, html, dash_table
import scanpy as sc
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from . import config
from .functions import *
import plotly.express as px
from plotly.subplots import make_subplots
from time import time
import json
import dash_cytoscape as cyto

def cytoscape_graph():
    return cyto.Cytoscape(
        id='graph_cytoscape',
        layout={'name': 'preset'},
        style={'width': '100%', 'height': '500px'},
        elements=config.graph,
        userZoomingEnabled=True,  # Disable zooming
        userPanningEnabled=False,  # Disable panning
        stylesheet=[
                    {
                        'selector': 'node',
                        'style': {
                            'label': 'data(id)',  # Show the node labels
                            'text-wrap': 'wrap',
                            # 'text-max-width': 30,
                            # 'text-justification': 'left',
                            'text-margin-x':-21,
                            'text-margin-y':-24,
                            'text-valign': 'bottom',  # Center the label vertically
                            'text-halign': 'right',  # Center the label horizontally
                            'background-color': 'white',
                            'opacity':'data(opacity)',
                            # 'shape': 'box',
                            'width': 24.5,
                            'height': 25,
                            'border-color': 'black',
                            'border-width': .0,
                            'shape': 'square',
                            'font-size':2,
                            'z-index-compare':'manual',
                            'z-index':0,
                            'background-fit': 'cover',
                            'background-image': 'data(image)'
                        },
                    },
                    {
                        'selector': 'edge',
                        'style': {
                            'curve-style':'straight',
                            'line-cap':'round',
                            'source-endpoint':'outside-to-node',
                            'source-arrow-shape': 'circle',
                            'target-arrow-shape': 'circle',
                            'mid-target-arrow-shape': 'triangle',
                            'width':1,
                            'arrow-scale': .35,
                            'z-index-compare':'manual',
                            'z-index':1,
                            'source-endpoint': '10.5 0', 
                            'target-endpoint': '-10.5 0',
                            # 'line-color':'black'
                        },
                    },
                    {
                        'selector': ':selected',
                        'style': {
                            'border-color': 'red'
                        }
                    },
                ],
    )

def save_graph():

    image_unselected(config.selected)
    file = f"{config.analysis_folder}/analysis.json"
    with open(file,"w") as outfile:
        json_object = json.dumps(config.graph, indent=4, cls=NpEncoder)
        outfile.write(json_object)
    image_selected(config.selected)

def save_adata():

    pos = get_node_pos(config.selected)
    file = f"{config.analysis_folder}/h5ad/{config.graph[pos]['data']['h5ad_file']}"
    config.adata.write(file)

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

def get_ancestors(name):
    l = []
    node = get_node(name)
    while node['data']['parameters']['input'] not in [None]:
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

def list_observables():
    l = list(config.adata.obs.columns.values)
    ancestors = get_nodes()
    for i in ancestors:
        if "obs" in i["data"].keys():
            for j in i["data"]["obs"].keys():
                l += [f"{i['data']['name']}---{j}"]

    return l

def get_observable(obs):

    l = list(config.adata.obs.columns.values)
    
    if obs in l:
        return config.adata.obs[obs].values
    
    ancestors = get_nodes()
    for i in ancestors:
        if "obs" in i["data"].keys():
            for j in i["data"]["obs"].keys():
                l = f"{i['data']['name']}---{j}"
                if obs == l:
                    return np.array(i["data"]["obs"][j])

    return None

def change_node_selected(name):
    image_unselected(config.selected)
    config.selected = name
    image_selected(config.selected)

def prevent_race(name,computed=True,method=True):

    node = get_node(config.selected)
    if not node['data']['computed'] and computed:
        raise PreventUpdate()

    if node['data']['method'] != name and method:
        raise PreventUpdate()
    
def is_computed():

    return get_node(config.selected)["data"]["computed"]