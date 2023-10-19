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
import dash_ag_grid as dag

from .graph import *

from app import app

ARGINPUT = {
            "input":"Dropdown",
            "name":"input",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options": {"function": "node_names(exclude_downstream_from_node=config.selected)"}
        }

ARGBATCH = {
            "input":"Dropdown",
            "name":"batch",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options": {"function":"[i for i,j in zip(config.adata.obs.columns.values, config.adata.obs.dtypes) if j in ['str','object','category','int']]"}
        }

from .methods.qc import *

def fvalue(value):
    val = {"val":value}
    if type(value) == dict:
        value = value.copy()
        if "function" in value.keys():
            m = f"val = {value['function']}"
            exec(m, globals(), val)

    return val["val"]

def fvisible(arg):

    if "visible" in arg.keys():
        return fvalue(arg["visible"])
    else:
        return True

def args_eval(args):
    d = {}

    for i in args:
        val = fvalue(i["value"])
        d[i["name"]] = val

    return d

def parameters_eval(args):
    
    if type(args) == list:
        for i,j in enumerate(args):
            args[i] = parameters_eval(j)
    elif type(args) == dict:
        if "function" in args.keys():
            return fvalue(args)
        else:
            for i,j in args.items():
                args[i] = parameters_eval(j)

    return args

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

def node_rm(name):

    l = []
    for i,node in enumerate(config.graph):
        if 'id' in node['data'].keys():
            if node['data']['id'] != name:

                if name == node['data']['parameters']['input']: #Remove input from cells that have this node as input
                    node['data']['parameters']['input'] = None
                    config.adata.uns[node['name']]['parameters']['input'] = None

                l.append(node)
        else:
            if node['data']['target'] != name and node['data']['source'] != name:
                l.append(node)

    #Remove info from adata
    del_adata_node(name)

    config.graph = l

    return

def set_output(args):

    #obs
    obs_key = f"{config.selected}--"
    if "obs" in args.keys():
        for i,j in args["obs"].items():
            i = f"{obs_key}{i}"
            config.adata.obs[i] = j

    #var
    var_key = f"{config.selected}--"
    if "var" in args.keys():
        for i,j in args["var"].items():
            i = f"{var_key}{i}"
            config.adata.var[i] = j

    #obsm
    if "obsm" in args.keys():
        config.adata.obsm[config.selected] = args["uns"]
            
    #uns
    if "uns" in args.keys():
        config.adata.uns[config.selected] = args["uns"]

def make_arguments(id, arglist, loaded_args={}, add_execution_button=True, add_header=True):

    arglist = parameters_eval(arglist)

    l = [
        dbc.Col(
            dbc.Row(html.H1("Additional arguments")),
        ),
        dbc.Col(
            dbc.Row(
                dbc.Button("Save Image",id="analysis_saveimage_button", class_name="btn btn-primary btn-sm"),
            )
        )
    ]
    if add_header:
        l = [
            dbc.Row(
                [
                    dbc.Col(
                        html.H1(config.selected, id="analysis_name")
                    ),
                    # dbc.Col(
                    #     dbc.Row(
                    #         dbc.Button("Rename",id="analysis_rename_button", class_name="btn btn-primary btn-sm"),
                    #     ),
                    #     width={"offset":7},
                    #     align="center"
                    # ),
                    dbc.Col(
                        dbc.Row(
                            dbc.Button("Delete",id="analysis_delete_button", style={"background-color":"red"}, class_name="btn btn-primary btn-sm"),
                        ),
                        width={"offset":7},
                        align="center"
                    ),
                ],
                align="center"
            )
        ]

    for i,arg in enumerate(arglist):

        arg = arg.copy()

        if fvisible(arg):

            for i,j in arg.items():
                arg[i] = fvalue(j)

            if loaded_args != {}:
                try:
                    value = loaded_args[arg["name"]]
                except:
                    value = None
                    loaded_args[arg["name"]] = None
            else:
                value = arg["value"]
        
            l.append(
                dbc.Tooltip(
                        arg["description"],
                        target=id+str(i),
                        placement="bottom",
                )
            )
            lab = html.Label(arg["name"],id=id+str(i),style={'text-align': 'right'})
            if arg["input"] == "Input":

                input = dbc.Input(id="analysis_"+str(arg["name"]),value=value,type=arg["type"])

            elif arg["input"] == "Dropdown":

                input = dcc.Dropdown(
                            id="analysis_"+str(arg["name"]),
                            options=arg["options"],
                            value=value,
                            # placeholder="Select a column",
                            clearable=arg["clearable"]
                        )
                
            elif arg["input"] == "BooleanSwitch":

                input = daq.BooleanSwitch(id="analysis_"+str(arg["name"]), on=value)

            # elif arg["input"] == "QCTable":

            #     input = [
            #         dbc.Row(
            #             html.Div(
            #                     id="analysis_"+str(arg["name"]),
            #                     children=qc_table(arg["value"]),
            #             ),
            #         ),
            #     ]

            elif arg["input"] == "AgTable":

                if "deleteRows" in arg.keys():
                    d = arg["deleteRows"]
                else:
                    d = False

                input = [
                    dbc.Row(
                        html.Div(
                                children= ag_table("analysis_"+str(arg["name"]), arg["header"], value, deleteRows=d),
                        ),
                    ),
                ]

                if "addRows" in arg.keys():

                    input += [
                        dbc.Row(
                            html.Div(
                                    children= dbc.Button(id="analysis_"+str(arg["name"])+"_button", children="Add"),
                            ),
                        ),
                    ]

            else:

                print("ERROR: NO METHOD")

            l.append(
                dbc.Row(
                    [
                        dbc.Col(
                            lab,
                            width=2
                            # width=text_width
                        ),
                        dbc.Col(
                            input,
                            # width=input_width
                        )
                    ],
                )
            )

    if add_execution_button:
        l.append(
                dbc.Button("Execute",id="analysis_execute_button")
        )

    return l

# def args2params(data):

#     params = {}
#     for i in data:
#         if i["input"] == "QCTable":
#             l = []
#             for j in i["value"]:
#                 if j["name"] != "":
#                     if "del" in j.keys():
#                         del j["del"]
#                     l.append(j)
#             params["name"] = l
#         else:
#             params["name"] = i["value"]

#     return {i['name']:i['value'] for i in data}

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


# def qc_table(rowData):

#     a = [" "]+list(config.adata.var.columns.values)

#     columnDefs = [
#         {
#             "headerName": "",
#             "field":"del", 
#             "cellRenderer": "Button", 
#             "lockPosition":'left', 
#             "cellRendererParams": {"className": "btn btn-danger"},  
#             "editable": False
#         },
#         {
#             "headerName": "Name metric",
#             "field": "name",
#             "editable": True
#         },
#         {
#             "headerName": ".var",
#             "field": "var",
#             "cellEditor": {"function": "DCC_Dropdown"},
#             "cellEditorParams": {"function": f"dynamicOptions(params,{a})"},
#             "cellEditorPopup": True,
#             "cellEditorPopupPosition": 'under',
#         },
#         {
#             "headerName": "Pattern",
#             "field": "pattern",
#             "editable": True
#         },
#         {
#             "headerName": "Style",
#             "field": "style",
#             "cellEditor": {"function": "DCC_Dropdown"},
#             "cellEditorParams": {"function": f"dynamicOptions(params,['counts','n_expressed_genes','proportion'])"},
#             "cellEditorPopup": True,
#             "cellEditorPopupPosition": 'under',
#         },
#         {
#             "headerName": "Genes",
#             "field": "genes",
#             "editable": False
#         },
#     ]

#     for i,d in enumerate(rowData):
#         if "del" not in d.keys():
#             rowData[i]["del"] = "Delete"
#         if "name" not in d.keys():
#             rowData[i]["name"] = ""
#         if "var" not in d.keys():
#             rowData[i]["del"] = a[0]
#         if "style" not in d.keys():
#             rowData[i]["pattern"] = "counts"
#         if "pattern" not in d.keys():
#             rowData[i]["pattern"] = ""
#         if "genes" not in d.keys():
#             rowData[i]["genes"] = "[]"

#     d = rowData[-1]
    
#     fig = dag.AgGrid(
#         id="analysis_measure",
#         rowData=rowData,
#         columnDefs=columnDefs,
#         defaultColDef={"editable": True},
#         deleteSelectedRows=True,
#         dashGridOptions={"suppressRowTransform":True, "suppressRowClickSelection":True},
#         columnSize="sizeToFit",
#     )

#     return fig

# @app.callback(
#     dash.Output("analysis_measure","rowData", allow_duplicate=True),
#     dash.Input("analysis_measure","cellValueChanged"),
#     dash.Input("analysis_measure","cellRendererData"),
#     dash.State("analysis_measure","rowData"),
#     prevent_initial_call = True
# )
# def f(b,click,c):

#     if b != None:
#         if b["data"]["var"] not in  [""," "] and b["data"]["pattern"] != "":
#             v =  [i for i in config.adata.var[b["data"]["var"]].values if i.startswith(b["data"]["pattern"])]
#             if len(v) != config.adata.shape[1]:
#                 c[b["rowIndex"]]["genes"] = v
#             else:
#                 c[b["rowIndex"]]["genes"] = "All"
#         else:
#             c[b["rowIndex"]]["genes"] = "All"

#     if click:
#         if click["colId"] == "del":
#             c.pop(click["rowIndex"])   

#     return c

def load_node(name):

    if name in config.adata.uns.keys():
        config.active_node_parameters = get_node(config.selected)["data"]["parameters"].copy()

    change_node_selected(name)

    if name == "Raw":

        l = []

    else:

        l = [
            dbc.Row(
                        id='graph_main',
                        children=layout(name),
                        style={"margin-bottom":"1cm"}        
            ),
        ]

    config.selected = name

    return l

def layout(name_analysis):

    # node_pars = get_node(name_analysis)['data']['parameters']
    method = get_node(name_analysis)['data']['method']
    args = get_node(name_analysis)['data']['parameters']
    l = make_arguments(method, deepcopy(config.methods[method]["args"]["execution"]), args)
    l2 = []
    p = []
    if get_node(name_analysis)['data']['computed']:
        l2 = make_arguments(method, deepcopy(config.methods[method]["args"]["postexecution"]), args, add_execution_button=False, add_header=False)
        p = config.methods[method]["plot"]()
        
    style = dict()
    style["background-color"]="lightgray"
    
    return  [
                dbc.Row(l,  
                        id = 'analysis_args',
                        # width='50%',
                        style=style
                        ),
                dbc.Row(id='analysis_postargs',
                        children=l2,  
                        # width='50%',
                        style=style
                        ),
                dbc.Row(id='analysis_plot',children=p)
            ]