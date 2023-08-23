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

from .methods.leiden import *
from .methods.louvain import *

from app import app

name_analysis = ""
functions = {}
names = ["Leiden","Louvain"]
methods = ["leiden","louvain"]

#Lists
def layout():
    return dbc.Container(
        [
            dbc.Row(
                children=[
                    dbc.Col(html.H1("Clustering"), width="auto")
                ],
                justify="center",
                className="mb-4"
            ),
            dbc.Row(
                children=[
                    dbc.Col(html.H1("Clusterings Performed"), width="auto")
                ],
                justify="left",
            ),
            dbc.Row(
                children=[
                dash_table.DataTable(
                    id='clustering_table_methods',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Name","Type"]
                    ],
                    data=make_clustering_table_methods_data(config.adata),
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
                            id='clustering_dropdown_methods',
                            value=None,
                            options=make_clustering_dropdown_methods(config.adata),
                            clearable=True,
                        )
                    ),
                    dbc.Col(
                        dbc.Button(
                            id='clustering_load_button',
                            children="Load"
                        )
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id='clustering_dropdown_analysis',
                            value=None,
                            options=names,
                            clearable=False,
                        )
                    ),
                    dbc.Col(
                        dbc.Button(
                            id='clustering_new_button',
                            children="New"
                        )
                    )
                ]
            ),
            dbc.Row(
                children=[
                    dbc.Col(html.H1(id='clustering_method', children=""), width="auto")
                ],
                justify="left",
            ),
            dbc.Row(
                id='clustering_analysis'
            ),
            dbc.Row(
                [
                    dbc.Button(id='clustering_save_button', children="Save",
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
            dbc.Row(
                children=[
                    dbc.Col(html.H1("Annotation"), width="auto")
                ],
                justify="center",
                className="mb-4"
            ),

        ]
    )

############################################################################################################################################
############################################################################################################################################
# Functions
############################################################################################################################################
############################################################################################################################################
def make_clustering_table_methods_data(adata):

    global methods

    if "__interactive__" in adata.uns.keys():
        return [{"Name":i,"Type":j["type"]} for i,j in adata.uns["__interactive__"].items() if j["type"] in methods]
    else:
        return []
    
def make_clustering_dropdown_methods(adata):

    global methods

    if "__interactive__" in adata.uns.keys():
        return [i for i,j in adata.uns["__interactive__"].items() if j["type"] in methods]
    else:
        return []
    
for name, method in zip(names,methods):

    callback_code = make_method_layout(method)
    exec(callback_code)

    add_function_to_dict = f"""functions["{name}"] = make_{method}"""
    exec(add_function_to_dict)

############################################################################################################################################
############################################################################################################################################
# Callbacks
############################################################################################################################################
############################################################################################################################################
@app.callback(
    dash.Output('clustering_method', 'children', allow_duplicate=True),
    dash.Output('clustering_analysis', 'children', allow_duplicate=True),
    [
        dash.Input('clustering_load_button', 'n_clicks')
    ],
    [
        dash.State('clustering_dropdown_methods', 'value')
    ],
    prevent_initial_call=True
)
def clustering_analysis(_, name):

    global name_analysis
    global functions
    name_analysis = name

    analysis = config.adata.uns["__interactive__"][name_analysis]["type"]

    l = [
        dbc.Row(
                    id='feature_selection_main',
                    children=functions[analysis](name_analysis, config.adata),
                    style={"margin-bottom":"1cm"}        
        ),
    ]

    return name, l

@app.callback(
    dash.Output('clustering_method', 'children', allow_duplicate=True),
    dash.Output('clustering_analysis', 'children', allow_duplicate=True),
    [
        dash.Input('clustering_new_button', 'n_clicks')
    ],
    [
        dash.State('clustering_dropdown_analysis', 'value')
    ],
    prevent_initial_call=True
)
def clustering_analysis(_, analysis):

    if "__interactive__" not in config.adata.uns.keys():
        config.adata.uns["__interactive__"] = {}

    name = analysis
    count = 0
    while name in config.adata.uns["__interactive__"].keys():
        count += 1
        name = analysis + " " + str(count)

    global name_analysis 
    global functions
    name_analysis = name
    
    config.adata.uns["__interactive__"][name_analysis] = {"type":analysis,"computed":False,"params":{}}

    l = [
        dbc.Row(
                    id='clustering_main',
                    children=functions[analysis](name_analysis, config.adata),
                    style={"margin-bottom":"1cm"}        
        ),
    ]

    return name, l

@app.callback(
    dash.Output('clustering_table_methods', 'data', allow_duplicate=True),
    dash.Output('clustering_dropdown_methods', 'options', allow_duplicate=True),
    [
        dash.Input('clustering_save_button', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def save_clustering(n_clicks):

    if config.CACHEFOLDER in config.file_path:
        config.adata.write(config.file_path)
    else:
        config.adata.write(config.CACHEFOLDER+config.file_path)

    return make_clustering_table_methods_data(config.adata), make_clustering_dropdown_methods(config.adata)

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
        add_input = "config.adata.uns['__interactive__'][name_analysis]['__input__'] = input"

    callback_code = f"""
@app.callback(
    dash.Output('{method}_plots1', 'children'),
    dash.Output('{method}_plots2', 'children'),
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

        config.adata.uns["__interactive__"][name_analysis]["computed"] = True
        {add_input}
        config.adata.uns["__interactive__"][name_analysis]["params"] = {{ {dictargs} }}

    return make_{method}_plots1(config.adata, name_analysis), make_{method}_plots2(config.adata, name_analysis)
"""

    exec(callback_code, globals(), locals())