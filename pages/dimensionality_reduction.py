import os
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import scanpy as sc
import dash_table
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from . import config
from .functions import *
import plotly.express as px
from plotly.subplots import make_subplots

from .methods.qc import *
from .methods.filtering import *
from .methods.feature_selection import *
from .methods.normalize import *
from .methods.log1p import *
from .methods.pca import *
from .methods.neighbors import *
from .methods.umap import *
from .methods.leiden import *
from .methods.louvain import *

import dash_cytoscape as cyto

from app import app

functions = {}
arguments = {}
names = ["QC","Filtering","Feature Selection","Normalize","Log1p","PCA","Neighbors","UMAP","Leiden","Louvain"]
methods = ["qc","filtering","feature_selection","normalize","log1p","pca","neighbors","umap","leiden","louvain"]

node_clicked = ""
old_load = 0

#Lists
def layout():
    return dbc.Container(
        [
            dbc.Row(
                children=[
                    dbc.Col(html.H1(id="title",children="Dimensionality Reduction"), width="auto")
                ],
                justify="center",
                className="mb-4"
            ),
            dbc.Row(
                children=[
                dash_table.DataTable(
                    id='dimensionality_reduction_table_methods',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Name","Type","Method"]
                    ],
                    data=make_dimensionality_reduction_table_methods_data(config.adata),
                    editable=False,
                    row_deletable=True,
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                ),
                ]
            ),
            dbc.Row(
                children=[
                    dbc.Col(
                        dcc.Dropdown(
                            id='dimensionality_reduction_dropdown_methods',
                            value=None,
                            options=make_dimensionality_reduction_dropdown_methods(config.adata),
                            clearable=True,
                        )
                    ),
                    dbc.Col(
                        dbc.Button(
                            id='dimensionality_reduction_load_button',
                            children="Load"
                        )
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id='dimensionality_reduction_dropdown_analysis',
                            value=names[0],
                            options=names,
                            clearable=False,
                        )
                    ),
                    dbc.Col(
                        dbc.Button(
                            id='dimensionality_reduction_new_button',
                            children="New"
                        )
                    )
                ]
            ),
            dbc.Row(
                cyto.Cytoscape(
                    id='dimensionality_reduction_cytoscape',
                    layout={'name': 'preset'},
                    style={'width': '100%', 'height': '500px'},
                    elements=make_cytoscape(config.adata),
                    userZoomingEnabled=True,  # Disable zooming
                    userPanningEnabled=False,  # Disable panning
                    stylesheet=[
                                {
                                    'selector': 'node',
                                    'style': {
                                        'label': 'data(summary)',  # Show the node labels
                                        'text-wrap': 'wrap',
                                        # 'text-max-width': 30,
                                        # 'text-justification': 'left',
                                        'text-margin-x':-19,
                                        'text-margin-y':-24,
                                        'text-valign': 'bottom',  # Center the label vertically
                                        'text-halign': 'right',  # Center the label horizontally
                                        'background-color': 'data(color)',
                                        # 'shape': 'box',
                                        'width': 20,
                                        'height': 25,
                                        'border-color': 'black',
                                        'border-width': .5,
                                        'shape': 'square',
                                        'font-size':2    
                                    },
                                },
                                {
                                    'selector': ':selected',
                                    'style': {
                                        'border-color': 'red'
                                    }
                                },
                            ],
                ),
                style={'width': '100%', 'margin': 'auto'},
            ),
            # ]),
            dbc.Modal(
                id='info-modal',
                size='lg',
                is_open=False,
                children=[
                    dbc.ModalHeader('Node Information'),
                    dbc.ModalBody(id='modal-body'),
                    dbc.ModalFooter(
                        [
                            dbc.Button('Load', id='load-modal', className='ml-auto', n_clicks=0),
                            dbc.Button('Close', id='close-modal', className='ml-auto', n_clicks=0)
                        ]
                    )
                ]
            ),
            dbc.Row(
                id='dimensionality_reduction_analysis'
            ),
            dbc.Row(
                [
                    dbc.Button(id='dimensionality_reduction_save_button', children="Save",
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
        ],
        fluid=True
    )

############################################################################################################################################
############################################################################################################################################
# Functions
############################################################################################################################################
############################################################################################################################################
def make_dimensionality_reduction_table_methods_data(adata):

    global methods

    if "__interactive__" in adata.uns.keys():
        return [{"Name":i,"Type":j["type"]} for i,j in adata.uns["__interactive__"].items() if j["type"] in names]
    else:
        return []

def make_dimensionality_reduction_dropdown_methods(adata):

    global methods

    if "__interactive__" in adata.uns.keys():
        return [i for i,j in adata.uns["__interactive__"].items() if j["type"] in names]
    else:
        return []

def make_cytoscape(adata,positions=[]):

    nodes = [
            {
                'data': {'id': 'Raw', 'label': 'Raw', 'color': 'white', 'summary':'Raw'}, 
                'position':{'x':-200,'y':0},
            }
    ]
    edges=[]

    colormap = {
        "QC":"lightgreen",
        "Filtering":"orange",
    }

    if "__interactive__" in adata.uns.keys():
        for i,j in adata.uns["__interactive__"].items():

            if j["computed"]:
                color = colormap[j["type"]]
            else:
                color = "white"

            if 'position' not in j.keys():
                j['position'] = {'x':0,'y':0}

            pos = j["position"]

            for node in positions:
                if node["data"]["id"] == i:
                    pos = node["position"]

            summary = f"{i}\nMethod:{j['method']}\n\n"
            for prop in arguments[j["method"]](adata):   
                if type(prop) != str:
                    if "summary" in prop.keys():
                        summary += f"{prop['name']}:\n "+str(j['params'][prop['name']].split('[')[-1].split(']')[0].replace(',','\n '))

            config.adata.uns["__interactive__"][i]["position"]=pos
            nodes.append(
                {
                    'data': {'id': str(i), 'label': f'{i}', 'color': color, 'type':j["type"] , 'method':j["type"], 'summary':summary, 'params':j['params']}, 
                    'position':pos,
                }
            )

        for i,j in adata.uns["__interactive__"].items():
            if "__input__" in j.keys():
                if j["__input__"] != None:
                    if j["__input__"] in adata.uns["__interactive__"].keys() or j["__input__"] == "Raw":
                        edges.append(
                            {'data': {'source': j["__input__"], 'target': str(i)}}
                        )

    return nodes+edges

for name, method in zip(names,methods):

    callback_code = make_method_layout(method)
    exec(callback_code)

    add_function_to_dict = f"""arguments["{name}"] = {method}_args"""
    exec(add_function_to_dict)

    add_function_to_dict = f"""functions["{name}"] = make_{method}"""
    exec(add_function_to_dict)

############################################################################################################################################
############################################################################################################################################
# Callbacks
############################################################################################################################################
############################################################################################################################################
@app.callback(
    dash.Output('dimensionality_reduction_analysis', 'children', allow_duplicate=True),
    [
        dash.Input('dimensionality_reduction_load_button', 'n_clicks')
    ],
    [
        dash.State('dimensionality_reduction_dropdown_methods', 'value')
    ],
    prevent_initial_call=True
)
def dimensionality_reduction_analysis(_, name):

    global functions
    config.name_analysis = name

    analysis = config.adata.uns["__interactive__"][config.name_analysis]["type"]

    l = [
        dbc.Row(
                    id='feature_selection_main',
                    children=functions[analysis](config.name_analysis, config.adata),
                    style={"margin-bottom":"1cm"}        
        ),
    ]

    return l

@app.callback(
    dash.Output('dimensionality_reduction_analysis', 'children', allow_duplicate=True),
    dash.Output('dimensionality_reduction_table_methods', 'data', allow_duplicate=True),
    dash.Output('dimensionality_reduction_dropdown_methods', 'options', allow_duplicate=True),
    [
        dash.Input('dimensionality_reduction_new_button', 'n_clicks')
    ],
    [
        dash.State('dimensionality_reduction_dropdown_analysis', 'value')
    ],
    prevent_initial_call=True
)
def dimensionality_reduction_analysis(_, analysis):

    if "__interactive__" not in config.adata.uns.keys():
        config.adata.uns["__interactive__"] = {}

    name = analysis
    count = 0
    while name in config.adata.uns["__interactive__"].keys():
        count += 1
        name = analysis + " " + str(count)
 
    global functions
    config.name_analysis = name
    
    params = {}
    for i in arguments[analysis](config.adata):
        if type(i) != str:
            params[i["name"]] = i["value"]

    config.adata.uns["__interactive__"][config.name_analysis] = \
        {
        "type":analysis,
        "method":analysis,
        "position":{"x":0,"y":0},
        "computed":False,
        "params":params
        }

    l = [
        dbc.Row(
                    id='dimensionality_reduction_main',
                    children=functions[analysis](config.name_analysis, config.adata),
                    style={"margin-bottom":"1cm"}        
        ),
    ]

    return l, make_dimensionality_reduction_table_methods_data(config.adata), make_dimensionality_reduction_dropdown_methods(config.adata)

@app.callback(
    dash.Output('dimensionality_reduction_dropdown_methods', 'value', allow_duplicate=True),
    dash.Input('dimensionality_reduction_cytoscape', 'tapNode'),
    prevent_initial_call=True
)
def display_click_data(tap_node_data):

    if tap_node_data is not None:
        global node_clicked
        node_clicked = tap_node_data['data']['label']
        return tap_node_data['data']['label']
    
    node_clicked = ""
    return None

@app.callback(
    dash.Output('dimensionality_reduction_table_methods', 'data', allow_duplicate=True),
    dash.Output('dimensionality_reduction_dropdown_methods', 'options', allow_duplicate=True),
    dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('modal-name','children'),
    dash.Input('modal-change_name_button','n_clicks'),
    dash.State('modal-change_name','value'),
    prevent_initial_call=True
)
def change_name(n_clicks, name):

    global node_clicked

    if name != node_clicked:
        config.adata.uns['__interactive__'][name] = config.adata.uns['__interactive__'][node_clicked].copy()
        del config.adata.uns['__interactive__'][node_clicked]
        node_clicked = name

    return make_dimensionality_reduction_table_methods_data(config.adata), \
            make_dimensionality_reduction_dropdown_methods(config.adata), \
            make_cytoscape(config.adata), \
            f"Name: {node_clicked}"

@app.callback(
    dash.Output('info-modal', 'is_open'),
    dash.Output('modal-body', 'children'),
    dash.Output('dimensionality_reduction_analysis', 'children', allow_duplicate=True),
    dash.Input('dimensionality_reduction_cytoscape', 'tapNode'),
    dash.Input('load-modal', 'n_clicks'),
    dash.Input('close-modal', 'n_clicks'),
    dash.State('info-modal', 'is_open'),
    dash.State('dimensionality_reduction_analysis', 'children'),
    dash.State('modal-body', 'children'),
    prevent_initial_call=True
)
def open_modal(tap_node, load_clicks, close_clicks, is_open, dra, name_dict):
    ctx = dash.callback_context
    triggered_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    global node_clicked

    if triggered_id == 'dimensionality_reduction_cytoscape' and tap_node:
        node_label = tap_node['data']['label']

        if node_label == "Raw":
            return False, None, dra
        
        l=[]
        for i,j in config.adata.uns["__interactive__"][node_label]["params"].items():
            l.append(dbc.Label(f"\t{i}={j}"))

        if node_clicked == node_label:
            modal_body = dbc.Row(
                [
                    dbc.Label(id='modal-name',children=f"Name: {node_label}"),
                    dbc.Col(
                        dbc.Input(id="modal-change_name",value=node_label),
                    ),
                    dbc.Col(
                        dbc.Button(id='modal-change_name_button',children="Change Name")
                    ),
                    dbc.Label(f"Type: {tap_node['data']['type']}"),
                    dbc.Label(f"Method: {tap_node['data']['method']}"),
                    dbc.Label("Arguments:"),
                ]+l
            )
            return True, modal_body, dra
        else:
            node_clicked = node_label
            return False, None, dra
    
    if triggered_id == 'close-modal' and close_clicks:
        return False, None, dra
    
    if triggered_id == 'load-modal' and load_clicks:

        name = name_dict['props']['children'][1]["props"]["children"]["props"]["value"]


        global functions
        config.name_analysis = name

        analysis = name_dict['props']['children'][4]['props']['children'].split("Method: ")[-1]

        l = [
            dbc.Row(
                        id='feature_selection_main',
                        children=functions[analysis](config.name_analysis, config.adata),
                        style={"margin-bottom":"1cm"}        
            ),
        ]

        return False, None, l
    
    return is_open, None

for method in methods:

    f = method+"_args"
    function_code = f"a = {f}(config.adata)"

    exec(function_code, globals(), locals())

    m = "["
    args = ""
    kwargs = ""
    dictargs = ""
    for i in a:
        if type(i) == dict:
            j = i["name"]
            if i["input"] == 'BooleanSwitch':
                m += f"dash.State('{method}_{j}', 'on'),"
            elif i["input"] == 'ThresholdTable':
                m += f"dash.State('{method}_{j}', 'data'),"
            elif i["input"] == 'MeasureTable':
                m += f"dash.State('{method}_{j}', 'children'),"
            else:
                m += f"dash.State('{method}_{j}', 'value'),"
            args += f"{j},"
            kwargs += f"{j}={j},"
            dictargs += f"'{j}':{j},"
    m = m[:-1]
    args = args[:-1]
    kwargs = kwargs[:-1]
    dictargs = dictargs[:-1]
    m += "]"
    add_input = ""
    if "input" in args:
        add_input = "config.adata.uns['__interactive__'][config.name_analysis]['__input__'] = input"

    callback_code = f"""
@app.callback(
    dash.Output('{method}_plots1', 'children', allow_duplicate=True),
    dash.Output('{method}_plots2', 'children', allow_duplicate=True),
    dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    [
        dash.Input('button_{method}', 'n_clicks'),
    ],
    dash.State('dimensionality_reduction_cytoscape', 'elements'),
    {m},
    prevent_initial_call=True
)
def function_{method}(_, pos_, {args}):
    if _ != None:

        f_{method}(config.adata, config.name_analysis, {kwargs})

        config.adata.uns["__interactive__"][config.name_analysis]["computed"] = True
        {add_input}
        config.adata.uns["__interactive__"][config.name_analysis]["params"] = {{ {dictargs} }}

    return make_{method}_plots1(config.adata, config.name_analysis), make_{method}_plots2(config.adata, config.name_analysis), make_cytoscape(config.adata, pos_)
"""

    exec(callback_code, globals(), locals())

    m = ""
    args = ""
    kwargs = ""
    for i in a:
        if type(i) == dict:
            j = i["name"]
            if i["input"] == 'BooleanSwitch':
                m += f"dash.Input('{method}_{j}', 'on'),"
            else:
                m += f"dash.Input('{method}_{j}', 'value'),"
            args += f"{j},"
            kwargs += f"{j}={j},"
    m = m[:-1]
    args = args[:-1]
    kwargs = kwargs[:-1]

    turn_off_code = f"""
@app.callback(
     dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    [
        dash.Input('button_{method}', 'n_clicks'),
        {m}
    ],
    dash.State('dimensionality_reduction_cytoscape', 'elements'),
    dash.State('dimensionality_reduction_load_button', 'n_clicks'),
    # dash.State('dimensionality_reduction_cytoscape', 'elements'),
    prevent_initial_call=True
)
def function_off_{method}(
    _,{args}, pos_, load_
):

    global old_load

    if config.adata.uns["__interactive__"][config.name_analysis]["params"] != dict():

        m = [i["name"] for i in {method}_args(config.adata) if type(i) == dict]
        l = True
        for i,j in zip(m,[{args}]):
            l2 = config.adata.uns["__interactive__"][config.name_analysis]["params"][i] == j
            ll = l and l2

        if not ll:
            print([config.adata.uns["__interactive__"][config.name_analysis]["params"][i] for i in m],{args})
            config.adata.uns["__interactive__"][config.name_analysis]["computed"] = False
            deactivate_downstream(config.adata, config.name_analysis)
        

    return make_cytoscape(config.adata, pos_)    
"""
    exec(turn_off_code, globals(), locals())

    m = f"dash.Input('{method}_input', 'value'),"

    link_code = f"""
@app.callback(
     dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    [
        {m}
    ],
    prevent_initial_call=True
)
def function_link_{method}(value):
    
    config.adata.uns["__interactive__"][config.name_analysis]["__input__"] = value

    return make_cytoscape(config.adata)
"""
    exec(link_code, globals(), locals())

@app.callback(
    dash.Output('dimensionality_reduction_table_methods', 'data', allow_duplicate=True),
    dash.Output('dimensionality_reduction_dropdown_methods', 'options', allow_duplicate=True),
    dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    [
        dash.Input('dimensionality_reduction_save_button', 'n_clicks'),
    ],
    dash.State('dimensionality_reduction_cytoscape', 'elements'),
    prevent_initial_call=True
)
def save_dimensionality_reduction(n_clicks,pos_):

    if config.CACHEFOLDER in config.file_path:
        config.adata.write(config.file_path)
    else:
        config.file_path = config.CACHEFOLDER+config.file_path
        config.adata.write(config.file_path)

    return make_dimensionality_reduction_table_methods_data(config.adata), make_dimensionality_reduction_dropdown_methods(config.adata), make_cytoscape(config.adata,pos_)

@app.callback(
    dash.Output('dimensionality_reduction_table_methods', 'data', allow_duplicate=True),
    dash.Output('dimensionality_reduction_dropdown_methods', 'options', allow_duplicate=True),
    dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('dimensionality_reduction_analysis', 'children', allow_duplicate=True),
    [
        dash.Input('dimensionality_reduction_table_methods', 'data'),
    ],
    [
        dash.State('dimensionality_reduction_analysis', 'children'),
        dash.State('dimensionality_reduction_cytoscape', 'elements'),
    ],
    prevent_initial_call=True
)
def remove_dimensionality_reduction(data, children, pos_):

    l = [i["Name"] for i in data]


    ls = [i for i,j in config.adata.uns["__interactive__"].items()]
    prop = [j for i,j in config.adata.uns["__interactive__"].items()]

    for i,prop in zip(ls,prop):
        j = config.adata.uns["__interactive__"][i]
        if (i not in l) and (prop["type"] in names):
            deactivate_downstream(config.adata, i)      
            for j,k in config.adata.uns["__interactive__"].items():
                if "__input__" in k.keys():
                    if k["__input__"] == i:
                        k["__input__"] = None
            try:
                del config.adata.uns["__interactive__"][i]
            except:
                None

            if config.name_analysis == i:
                children = []

            if config.name_analysis == i:
                children2 = ""


    if config.adata.uns["__interactive__"] == {}:
        del config.adata.uns["__interactive__"]

    return make_dimensionality_reduction_table_methods_data(config.adata), make_dimensionality_reduction_dropdown_methods(config.adata), make_cytoscape(config.adata,pos_), children84