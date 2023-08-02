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

from .args.feature_selection import *
from .args.pca import *
from .args.umap_args import umap_args
from .args.neighbors_args import neighbors_args
import dash_cytoscape as cyto

from app import app

name_analysis = ""
functions = {}

#Lists
def layout():
    return dbc.Container(
        [
            dbc.Row(
                children=[
                    dbc.Col(html.H1("Dimensionality Reduction"), width="auto")
                ],
                justify="center",
                className="mb-4"
            ),
            dbc.Row(
                children=[
                    dbc.Col(html.H1("Dimensionality Reduction Performed"), width="auto")
                ],
                justify="left",
            ),
            dbc.Row(
                children=[
                dash_table.DataTable(
                    id='dimensionality_reduction_table_methods',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Name","Type"]
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
                            value=None,
                            options=['Feature Selection','PCA','Neighbors','UMAP'],
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
                    style={'width': '80%', 'height': '500px'},
                    elements=make_cytoscape(config.adata),
                    userZoomingEnabled=False,  # Disable zooming
                    userPanningEnabled=False,  # Disable panning
                    stylesheet=[
                                {
                                    'selector': 'node',
                                    'style': {
                                        'label': 'data(label)',  # Show the node labels
                                        'text-valign': 'top',  # Center the label vertically
                                        'text-halign': 'center',  # Center the label horizontally
                                        'background-color': 'data(color)',
                                        'shape': 'circle',
                                        'width': 50,
                                        'height': 50,
                                        'border-color': 'black',
                                        'border-width': 1
                                    }
                                },
                                {
                                    'selector': ':selected',
                                    'style': {
                                        'border-color': 'red'
                                    }
                                }
                            ]
                ),
                style={'width': '90%', 'margin': 'auto'},
            ),
            dbc.Row(
                children=[
                    dbc.Col(html.H1(id='dimensionality_reduction_method', children="Method: "), width="auto")
                ],
                justify="left",
            ),
            dbc.Row(
                id='dimensionality_reduction_analysis'
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
    if "dimensionality_reduction" in adata.uns.keys():
        return [{"Name":i,"Type":j["type"]} for i,j in adata.uns["dimensionality_reduction"].items()]
    else:
        return []

def make_dimensionality_reduction_dropdown_methods(adata):
    if "dimensionality_reduction" in adata.uns.keys():
        return [i for i in adata.uns["dimensionality_reduction"].keys()]
    else:
        return []

def make_cytoscape(adata):

    xmap = {
            "Feature Selection":0,
            "PCA":1,
            "Neighbors":2,
            "UMAP":3            
            }
    ymap = {
            "Feature Selection":0,
            "PCA":0,
            "Neighbors":0,
            "UMAP":0            
    }

    if "dimensionality_reduction" in adata.uns.keys():
        nodes = []
        for i,j in adata.uns["dimensionality_reduction"].items():
            if j["__computed__"]:
                color = 'green'
            else:
                color = 'lightgray'

            nodes.append(
                {
                    'data': {'id': str(i), 'label': f'{i}', 'color': color}, 
                    'position':{'x':xmap[j["type"]]*200,'y':ymap[j["type"]]*100},
                }
            )
            ymap[j["type"]] += 1    

        edges=[]
        for i,j in adata.uns["dimensionality_reduction"].items():
            if "input" in j.keys():
                if j["input"] != None:
                    edges.append(
                        {'data': {'source': j["input"], 'target': str(i)}}
                    )

        return nodes+edges
    else:
        return []

for name, method in zip(["Feature Selection","PCA","Neighbors"],["feature_selection","pca","neighbors"]):

    callback_code = f"""
def make_{method}(name_analysis, adata):

    l = make_arguments('{method}', {method}_args(adata), adata.uns["dimensionality_reduction"][name_analysis]["params"])
        
    style = dict()
    style["background-color"]="lightgray"
    
    return  [
                dbc.Col(l,
                        width=4,
                        style=style
                        ),
                dbc.Col([
                    dbc.Row(id='{method}_plots1',children=make_{method}_plots1(adata, name_analysis))
                ]),
                dbc.Row(id='{method}_plots2',children=make_{method}_plots2(adata, name_analysis))
            ]

functions["{name}"] = make_{method}
"""

    exec(callback_code)

def make_neighbors_plots1(adata):
    return []

def make_neighbors_plots2(adata):
    return []

def make_umap(name_analysis, adata):

    l = make_arguments("umap",umap_args(adata))

    return  [
                dbc.Col(l,width=4,
                        style={"background-color":"lightgray"}
                        ),
                dbc.Col([
                    dbc.Row(id="umap_plot",children=make_umap_plots(adata))
                ]
                ),
            ]

functions["UMAP"] = make_umap

def make_umap_plots(adata):

    global name_analysis

    return []

    # if "X_"+name_analysis in adata.obsm.keys():

    #     x = adata.obsm[:,0]
    #     y = adata.obsm[:,1]

    #     return [
    #         dbc.Col(
    #             dcc.Graph(
    #                 figure = {'data':[
    #                             go.Scatter(
    #                                 x=x,
    #                                 y=y,
    #                                 mode='markers',
    #                                 name='Min threshold',
    #                             )]
    #                         }
    #             )
    #         ),
    #     ]
    # else:
    #     return []

############################################################################################################################################
############################################################################################################################################
# Callbacks
############################################################################################################################################
############################################################################################################################################
@app.callback(
    dash.Output('dimensionality_reduction_method', 'children', allow_duplicate=True),
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

    global name_analysis
    global functions
    name_analysis = name

    analysis = config.adata.uns["dimensionality_reduction"][name_analysis]["type"]

    l = [
        dbc.Row(
                    id='feature_selection_main',
                    children=functions[analysis](name_analysis, config.adata),
                    style={"margin-bottom":"1cm"}        
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
    ]

    return name, l

@app.callback(
    dash.Output('dimensionality_reduction_method', 'children', allow_duplicate=True),
    dash.Output('dimensionality_reduction_analysis', 'children', allow_duplicate=True),
    [
        dash.Input('dimensionality_reduction_new_button', 'n_clicks')
    ],
    [
        dash.State('dimensionality_reduction_dropdown_analysis', 'value')
    ],
    prevent_initial_call=True
)
def dimensionality_reduction_analysis(_, analysis):

    if "dimensionality_reduction" not in config.adata.uns.keys():
        config.adata.uns["dimensionality_reduction"] = {}

    name = analysis
    count = 0
    while name in config.adata.uns["dimensionality_reduction"].keys():
        count += 1
        name = analysis + " " + str(count)

    global name_analysis 
    global functions
    name_analysis = name
    
    config.adata.uns["dimensionality_reduction"][name_analysis] = {"type":analysis,"__computed__":False,"params":{}}

    l = [
        dbc.Row(
                    id='feature_selection_main',
                    children=functions[analysis](name_analysis, config.adata),
                    style={"margin-bottom":"1cm"}        
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
    ]

    return name, l

@app.callback(
    dash.Output('dimensionality_reduction_dropdown_methods', 'value', allow_duplicate=True),
    dash.Input('dimensionality_reduction_cytoscape', 'tapNode'),
    prevent_initial_call=True
)
def display_click_data(tap_node_data):

    if tap_node_data is not None:
        return tap_node_data['data']['label']
    
    return None

for method in ["feature_selection","pca","neighbors"]:

    f = method+"_args"
    function_code = f"a = {f}(config.adata)"

    exec(function_code, globals(), locals())

    m = "["
    args = ""
    kwargs = ""
    for i in a:
        if type(i) == dict:
            j = i["name"]
            if i["input"] == 'BooleanSwitch':
                m += f"dash.State('{method}_{j}', 'on'),"
            else:
                m += f"dash.State('{method}_{j}', 'value'),"
            args += f"{j},"
            kwargs += f"{j}={j},"
    m = m[:-1]
    args = args[:-1]
    kwargs = kwargs[:-1]
    m += "]"

    callback_code = f"""
@app.callback(
    dash.Output('{method}_plots1', 'children'),
    dash.Output('{method}_plots2', 'children'),
    dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    [
        dash.Input('button_{method}', 'n_clicks'),
    ],
    {m},
    prevent_initial_call=True
)
def function_{method}(_, {args}):
    if _ != None:
    
        global name_analysis

        f_{method}(config.adata, name_analysis, {kwargs})

        config.adata.uns["dimensionality_reduction"][name_analysis]["__computed__"] = True


    return make_{method}_plots1(config.adata, name_analysis), make_{method}_plots2(config.adata, name_analysis), make_cytoscape(config.adata)
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
    prevent_initial_call=True
)
def function_off_{method}(
    _,{args}
):
    
    global name_analysis

    if "params" in config.adata.uns["dimensionality_reduction"][name_analysis].keys():   

        if config.adata.uns["dimensionality_reduction"][name_analysis]["params"] != dict():

            m = [i["name"] for i in {method}_args(config.adata) if type(i) == dict]
            l = True
            for i,j in zip(m,[{args}]):
                l2 = config.adata.uns["dimensionality_reduction"][name_analysis]["params"][i] == j
                l = l and l2

            if not l:
                config.adata.uns["dimensionality_reduction"][name_analysis]["__computed__"] = False

    return make_cytoscape(config.adata)
"""
    exec(turn_off_code, globals(), locals())

@app.callback(
    dash.Output('dimensionality_reduction_table_methods', 'data', allow_duplicate=True),
    dash.Output('dimensionality_reduction_dropdown_methods', 'options', allow_duplicate=True),
    dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    [
        dash.Input('dimensionality_reduction_save_button', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def save_dimensionality_reduction(n_clicks):

    config.adata.write(config.file_path.split(".h5ad")[0]+"_qc.h5ad")

    return make_dimensionality_reduction_table_methods_data(config.adata), make_dimensionality_reduction_dropdown_methods(config.adata), make_cytoscape(config.adata)

@app.callback(
    dash.Output('dimensionality_reduction_table_methods', 'data', allow_duplicate=True),
    dash.Output('dimensionality_reduction_dropdown_methods', 'options', allow_duplicate=True),
    dash.Output('dimensionality_reduction_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('dimensionality_reduction_analysis', 'children', allow_duplicate=True),
    dash.Output('dimensionality_reduction_method', 'children', allow_duplicate=True),
    [
        dash.Input('dimensionality_reduction_table_methods', 'data'),
    ],
    [
        dash.State('dimensionality_reduction_analysis', 'children'),
        dash.State('dimensionality_reduction_method', 'children'),
    ],
    prevent_initial_call=True
)
def remove_dimensionality_reduction(data, children, children2):

    l = [i["Name"] for i in data]

    global name_analysis

    ls = list(config.adata.uns["dimensionality_reduction"].keys())

    for i in ls:
        j = config.adata.uns["dimensionality_reduction"][i]
        if i not in l:
            try:
                del config.adata.uns["dimensionality_reduction"][i]
            except:
                None

            if j["type"] == "Feature Selection":

                try:
                    config.adata.var.drop(i+"_highly_variable",inplace=True)
                except:
                    None
                    
            elif j["type"] == "PCA":

                try:
                    del config.adata.obsm["X_"+i]
                except:
                    None

            if name_analysis == i:
                children = []

            if name_analysis == i:
                children2 = ""


    if config.adata.uns["dimensionality_reduction"] == {}:
        del config.adata.uns["dimensionality_reduction"]

    return make_dimensionality_reduction_table_methods_data(config.adata), make_dimensionality_reduction_dropdown_methods(config.adata), make_cytoscape(config.adata), children, children2