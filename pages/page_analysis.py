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

from app import app

from .graph import *
from .arguments import *

from .methods.qc import *
from .methods.scrublet import *
from .methods.filtering import *
from .methods.log1p import *
from .methods.normalize import *
from .methods.feature_selection import *
from .methods.pca import *
from .methods.neighbors import *
from .methods.umap import *
from .methods.leiden import *
from .methods.louvain import *
from .methods.differential_expression import *
from .methods.marker_genes import *

methods = {
    "qc":{"method":"qc","type":"QC","recompute":False},
    # "scrublet":{"method":"scrublet","type":"QC","recompute":False},
    # "filtering":{"method":"filtering","type":"QC","recompute":True},
    # "log1p":{"method":"log1p","type":"Transformations","recompute":True},
    # "normalize":{"method":"normalize","type":"Transformations","recompute":True},
    # "feature_selection":{"method":"feature_selection","type":"DR","recompute":False},
    # "pca":{"method":"pca","type":"DR","recompute":False},
    # "neighbors":{"method":"neighbors","type":"Graph","recompute":False},
    # "umap":{"method":"umap","type":"Visualization","recompute":False},
    # "leiden":{"method":"leiden","type":"Clustering","recompute":False},
    # "louvain":{"method":"louvain","type":"Clustering","recompute":False},
    # "differential_expression":{"method":"differential_expression","type":"Clustering","recompute":False},
    # "marker_genes":{"method":"marker_genes","type":"Clustering","recompute":False},
}

graph_colormap={
    'QC':'yellow',
    # 'Transformations':'yellow',
    # 'DR':'yellow',
    # 'Graph':'yellow',
    # 'Visualization':'yellow',
    # 'Clustering':'yellow',
}

#Lists
def layout():
    return dbc.Container(
        [
            dbc.Row(
                children=[
                    dbc.Col(html.H1(id="title",children="Graphical Analysis"), width="auto"),
                    dbc.Col(),
                    dbc.Col(),
                    dbc.Col(
                        dbc.Button(id='graph_save_button', children="Save",
                                    size="lg",
                                    style={
                                        "background-color":"#343A40",
                                        'width': '280px', 
                                    }      
                                    )
                    )
                ],
                justify="center",
                className="mb-4"
            ),
            html.P(id="dumb",children=""),
            dbc.Modal(
                [
                    dbc.ModalHeader("Warning",close_button=False),
                    dbc.ModalBody(id="delete-message",
                                children=[]
                    ),
                    dbc.ModalFooter([
                            dbc.Button("Delete", id="delete-proceed", className="ml-auto"),
                            dbc.Button("Cancel", id="delete-cancel", className="ml-auto")
                    ])
                ],
                backdrop=False,
                id="delete-modal",
                size="sm",
            ),
            dbc.Modal(
                [
                    dbc.ModalHeader("Rename analysis",close_button=False),
                    dbc.ModalBody(id="rename-message",
                                children=[]
                    ),
                    dbc.ModalFooter([
                            dbc.Button("Rename", id="rename-proceed", className="ml-auto"),
                            dbc.Button("Cancel", id="rename-cancel", className="ml-auto")
                    ])
                ],
                backdrop=False,
                id="rename-modal",
                size="sm",
            ),
            dbc.Modal(
                [
                    dbc.ModalHeader("Warning",close_button=False),
                    dbc.ModalBody(children=html.Div("An input must be assigned before executing the cell.")),
                ],
                backdrop=True,
                id="execute-input-modal",
                size="sm",
            ),
            dbc.Modal(
                [
                    dbc.ModalHeader("Warning",close_button=False),
                    dbc.ModalBody(children=html.Div("The previous cell must be computed before computing this one.")),
                ],
                backdrop=True,
                id="execute-computed-modal",
                size="sm",
            ),
            dbc.Row(
                children=[
                dash_table.DataTable(
                    id='graph_table',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Name","Type","Method"]
                    ],
                    data=graph2table(),
                    editable=False,
                    row_deletable=False,
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                ),
                ]
            ),
            dbc.Row(
                children=[
                    dbc.Col(
                        dcc.Dropdown(
                            id='graph_dropdown_load',
                            value=None,
                            options=node_names(),
                            clearable=True,
                        )
                    ),
                    dbc.Col(
                        dbc.Button(
                            id='graph_load_button',
                            children="Load"
                        )
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id='graph_dropdown_analysis',
                            value=get_method("method")[0],
                            options=get_method("method"),
                            clearable=False,
                        )
                    ),
                    dbc.Col(
                        dbc.Button(
                            id='graph_new_button',
                            children="New"
                        )
                    )
                ]
            ),
            dbc.Row(
                cytoscape_graph(),
                style={'width': '100wh', 'height':'120wh', 'margin': 'auto'}, #vh/wh
            ),
            # ]),
            html.Div(id='graph_analysis',children=[]),
        ],
        fluid=True
    )

def get_method(name):
    global methods
    if name != "name":
        return [i[name] for _,i in methods.items()]
    else:
        return list(methods.keys())

def graph2table():
    return [{"Name":i['data']['id'],"Type":i['data']['type'],"Method":i['data']['method']} for i in get_nodes() if i['data']['id'] != 'Raw']

def args2summary(data):

    summary = ""
    for i in data:
        if 'summary' in i.keys():
            summary += f"{i['name']}:{str(i['value'])}"

    return summary

def make_interactive():
    if 'sc_interactive' not in config.adata.uns.keys():
        config.adata.uns['sc_interactive'] = {}

#Add new parameter
methods_implemented = []
for method in config.methods.keys():

    for i in config.methods[method]["args"]["execution"].copy():

        m_i = (i['name'],i['input'])

        if m_i not in methods_implemented:

            methods_implemented.append(m_i)

            input = f"dash.Input('analysis_{i['name']}','value'),"
            if i['input'] == 'BooleanSwitch':
                input = f"dash.Input('analysis_{i['name']}','on'),"
            elif i['input'] == 'QCTable':
                input = f"dash.Input('analysis_{i['name']}','rowData'),"

            add_function = f"""
@app.callback(
    dash.Output("dumb","children", allow_duplicate=True),
    {input}
    prevent_initial_call=True
)
def change_parameter_{i['name']}(data):

    config.active_node_parameters['{i['name']}'] = data

    return ""
"""

            exec(add_function, globals(), locals())

methods_implemented = []
for method in config.methods.keys():

    for i in config.methods[method]["args"]["postexecution"].copy():

        m_i = (i['name'],i['input'])

        if m_i not in methods_implemented:

            methods_implemented.append(m_i)

            input = f"dash.Input('analysis_{i['name']}','value'),"
            if i['input'] == 'BooleanSwitch':
                input = f"dash.Input('analysis_{i['name']}','on'),"
            elif i['input'] == 'QCTable':
                input = f"dash.Input('analysis_{i['name']}','rowData'),"

            add_function = f"""
@app.callback(
    dash.Output("analysis_postargs","children", allow_duplicate=True),
    dash.Output("analysis_plot","children", allow_duplicate=True),
    {input}
    prevent_initial_call=True
)
def change_parameter_{i['name']}(data):

    if data == None:
        raise PreventUpdate()

    config.adata.uns[config.selected]['parameters']['{i['name']}'] = data

    method = get_node(config.selected)["data"]["method"]
    post_args = config.methods[method]["args"]["postexecution"].copy()
    post_args_object = make_arguments(method, post_args, loaded_args=config.adata.uns[config.selected]['parameters'], add_execution_button=False, add_header=False)

    plot = config.methods[method]['plot']()

    return post_args_object, plot
"""

            exec(add_function, globals(), locals())


############################################################################################################################################
############################################################################################################################################
# Callbacks
############################################################################################################################################
############################################################################################################################################

#Update Dropdown (change of cytoscape)
@app.callback(
    dash.Output('graph_dropdown_load', 'options', allow_duplicate=True),
    [
        dash.Input('graph_cytoscape', 'elements')
    ],
    prevent_initial_call=True
)
def graph_update_dropdown(_):

    return node_names()

#Update Table (change of cytoscape)
@app.callback(
    dash.Output('graph_table', 'data', allow_duplicate=True),
    [
        dash.Input('graph_cytoscape', 'elements')
    ],
    prevent_initial_call=True
)
def graph_update_table(_):

    return graph2table()

@app.callback(
    dash.Output('dumb', 'children', allow_duplicate=True),
    dash.Input('graph_cytoscape', 'tapNode'),
    prevent_initial_call=True
)
def graph_update_pos(element):

    node_update_pos(config.graph, element['data']['id'], element['position'])
    config.max_x = np.max([i['position']['x'] for i in get_nodes()])

    return ""

#Add node
@app.callback(
    dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('graph_table', 'data', allow_duplicate=True),
    dash.Output('graph_analysis', 'children', allow_duplicate=True),
    [
        dash.Input('graph_new_button', 'n_clicks')
    ],
    [
        dash.State('graph_dropdown_analysis', 'value'),
        dash.State('graph_cytoscape', 'layout'),
    ],
    prevent_initial_call=True
)
def graph_new_node(_, method, cytoscape):

    #Make new name

    name = method
    nodes = node_names()
    count = 0
    while name in nodes:
        count += 1
        name = method + " " + str(count)

    args = config.methods[method]["args"]["execution"].copy()

    node = {
        'data': {
            'id': name, 
            'name': method,
            'method':methods[method]['method'], 
            'type':methods[method]['type'], 
            'color': graph_colormap[methods[method]['type']], 
            'image': '',
            'computed':False,
            'opacity':.3,
            'summary':'', 
            'parameters':args2params(args)
        }, 
        'position':{'x':config.max_x + 30,'y':0},
    }
    #make active node that is the one that will be presented
    # nodes_update_pos(config.graph,cytoscape)

    config.active_node_parameters = args2params(args).copy()
    config.graph.append(node)

    l = load_node(name)
    # make_node_summary(config.selected)

    return config.graph, graph2table(), l

# #Show edge when activate
# @app.callback(
#     dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
#     [
#         dash.Input('analysis_input', 'value')
#     ],
#     prevent_initial_call=True
# )
# def graph_edge(value):

#     for i,val in enumerate(config.graph):
#         #Change node input
#         if 'id' in val['data'].keys(): #Check is a node
#             if val['data']['id'] == config.selected:
#                 if val['data']['parameters']['input'] != value:
#                     config.graph[i]['data']['computed'] = False
#                     config.graph[i]['data']['opacity'] = .3
#                     config.graph[i]['data']['parameters']['input'] = value

#     for i,val in enumerate(config.graph):
#         #Change edge input
#         if 'target' in val['data'].keys(): #Check is an edge
#             if val['data']['target'] == config.selected: #there only one output to each node
#                 val['data']['source'] = value

#                 return config.graph

#     if value != None:
#         config.graph.append(
#             {
#                 'data':{'source':value,'target':config.selected},
#              }
#             )

#     return config.graph

#Execute analysis button
@app.callback(
    dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('execute-input-modal', 'is_open', allow_duplicate=True),
    dash.Output('execute-computed-modal', 'is_open', allow_duplicate=True),
    dash.Output('analysis_postargs', 'children', allow_duplicate=True),
    dash.Output('analysis_plot', 'children', allow_duplicate=True),
    [
        dash.Input('analysis_execute_button', 'n_clicks')
    ],
    dash.State('execute-input-modal', 'is_open'),
    dash.State('execute-computed-modal', 'is_open'),
    dash.State('analysis_plot', 'children'),
    prevent_initial_call=True
)
def execute(n_clicks, warning_input, warning_computed, plot):

    if n_clicks != None:

        if "input" in config.active_node_parameters.keys():

            if config.active_node_parameters['input'] == None:

                warning_input = True

            elif not get_node(config.active_node_parameters['input'])['data']['computed']:

                warning_computed = True

            else:

                #Remove previous edge
                make_interactive()
                node = get_node(config.selected)
                edge_rm(node['data']['parameters']['input'], config.selected)

                #Update node info
                update_node(config.selected)
                # make_node_summary(config.selected)

                node = get_node(config.selected)
                pos = get_node_pos(config.selected)
                params = get_node_parameters(config.selected, str2list=True)

                edge_add(params['input'], config.selected)

                inputArgs = get_args(params["input"], config.adata)

                outputArgs = config.methods[node['data']['method']]["function"](config.adata, inputArgs, params)

                if "uns" not in outputArgs.keys():
                    outputArgs["uns"] = {}
                outputArgs["uns"]["parameters"] = params

                set_output(config.adata, config.selected, outputArgs)

                #Add post arguments
                args = get_args(config.selected, config.adata)
                post_args = config.methods[method]["args"]["postexecution"].copy()
                for i in post_args:
                    val = fvalue(i["value"])
                    config.adata.uns[config.selected]["parameters"][i["name"]] = val

                activate_node(config.selected)
                deactivate_downstream(config.selected)

    post_args_object = []
    plot = []
    if get_node(config.selected)['data']['computed']:
        post_args = config.methods[method]["args"]["postexecution"].copy()
        post_args_object = make_arguments(method, post_args, add_execution_button=False, add_header=False)

        args = get_args(config.selected, config.adata)

        # print(args)
        plot = config.methods[node['data']['method']]["plot"]()

    return config.graph, warning_input, warning_computed, post_args_object, plot

#Delete analysis button
@app.callback(
    dash.Output('delete-modal', 'is_open'),
    dash.Output('delete-message', 'children'),
    dash.Input('analysis_delete_button', 'n_clicks'),
    prevent_initial_call=True
)
def delete(n_clicks):
        
    if n_clicks != None:

        l = [
                html.Div("You are about to remove an analysis node. This will remove all the computations performed over this node and descending nodes of the analysis. Are you sure that you want to proceed?"),
                html.Div(),
                html.Div(f"Analysis to be removed: {config.selected}")
            ]

        return True, l
    
    return False, []

#Proceed modal delete
@app.callback(
    dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('graph_analysis', 'children', allow_duplicate=True),
    dash.Output('delete-modal', 'is_open', allow_duplicate=True),
    [
     dash.Input('delete-proceed', 'n_clicks'),
    ],
    dash.State('graph_analysis', 'children'),
    prevent_initial_call=True
)
def delete_confirmation(n_clicks, analysis):
    
    modal = True
    if n_clicks != None:

        name = config.selected
        deactivate_downstream(name)
        node_rm(name)

        config.selected = 'Raw'
        analysis = ""
        modal = False
    
    return config.graph, analysis, modal

#Cancel modal delete
@app.callback(
    dash.Output('delete-modal', 'is_open', allow_duplicate=True),
    [
     dash.Input('delete-cancel', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def delete_cancel(n_clicks):
    
    modal = True
    if n_clicks != None:

        modal = False
    
    return modal

#Rename analysis button
@app.callback(
    dash.Output('rename-modal', 'is_open'),
    dash.Output('rename-message', 'children', allow_duplicate=True),
    dash.Input('analysis_rename_button', 'n_clicks'),
    prevent_initial_call=True
)
def rename(n_clicks):
        
    if n_clicks != None:

        l = [
            html.Div(f"Current name of analysis:\n"),
            html.Div(f"\t{config.selected}"),
            html.Div(f"Change by:"),
            dbc.Input(id="rename_name",value=config.selected,type="text")
        ]

        return True, l
    
    return False, []

#Proceed modal rename
@app.callback(
    dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('rename-modal', 'is_open', allow_duplicate=True),
    dash.Output('analysis_name', 'children', allow_duplicate=True),
    dash.Output('rename-message', 'children', allow_duplicate=True),
    [
     dash.Input('rename-proceed', 'n_clicks'),
    ],
    dash.State('rename_name', 'value'),
    dash.State('rename-message', 'children'),
    prevent_initial_call=True
)
def rename_confirmation(n_clicks, name, l):
    
    modal = True
    if n_clicks != None:

        n = node_names()

        if name not in n:

            node = get_node(config.selected)
            config.functions_method_rename[node['data']['method']](config.selected, name)
            node_rename(config.selected, name) #Rename graph
            config.selected = name #Rename configuration

            modal = False        
        else:
            l = [
                html.Div(f"THERE IS ALREADY AN ANALYSIS NODE WITH THIS NAME. CHOOSE OTHER NAME.\n", style={"background-color":"orange"}),
                html.Div(f"Current name of analysis:\n"),
                html.Div(f"\t{config.selected}"),
                html.Div(f"Change by:"),
                dbc.Input(id="rename_name",value=config.selected,type="text")
            ]    

    return config.graph, modal, config.selected, l

#Cancel modal rename
@app.callback(
    dash.Output('rename-modal', 'is_open', allow_duplicate=True),
    [
     dash.Input('rename-cancel', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def rename_cancel(n_clicks):
    
    modal = True
    if n_clicks != None:

        modal = False
    
    return modal

#Load button
@app.callback(
    dash.Output('graph_analysis', 'children', allow_duplicate=True),
    [
        dash.Input('graph_load_button', 'n_clicks')
    ],
    [
        dash.State('graph_dropdown_load', 'value'),
        dash.State('graph_analysis', 'children')
    ],
    prevent_initial_call=True
)
def graph_analysis(_, name, l):

    if _ != None and name != None:

        l = load_node(name)

    return l

#Double click to load node
@app.callback(
    dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('graph_analysis', 'children', allow_duplicate=True),
    dash.Input('graph_cytoscape', 'tapNode'),
    dash.State('graph_analysis', 'children'),
    prevent_initial_call=True
)
def display_click_data(tap_node_data, l):

    if tap_node_data is not None:

        t = time()
        name = tap_node_data['data']['id']

        if abs(config.tclick-t) < 1: #Double click

            l = load_node(name)
        
        else: # If not load, create again the buttons to avoid spurious reloadings

            if len(l) > 0: #Otherwise is raw
                l[0]["props"]["children"][0]["props"]["children"][0]["props"]["children"][1] = dbc.Col(
                        dbc.Row(
                            dbc.Button("Rename",id="analysis_rename_button", class_name="btn btn-primary btn-sm"),
                        ),
                        width={"offset":7},
                        align="center"
                    )
                l[0]["props"]["children"][0]["props"]["children"][0]["props"]["children"][2] = dbc.Col(
                        dbc.Row(
                            dbc.Button("Delete",id="analysis_delete_button", style={"background-color":"red"}, class_name="btn btn-primary btn-sm"),
                        ),
                        width={"offset":1},
                        align="center"
                    )
                l[0]["props"]["children"][-2]["props"]["children"][-1] = dbc.Button("Execute",id="analysis_execute_button")

            else:
                l = []

        config.tclick = t

        return config.graph, l

    return config.graph, l

#Save analysis
@app.callback(
    dash.Output("dumb","children",allow_duplicate=True),
    dash.Input("graph_save_button","n_clicks"),
    prevent_initial_call=True
)
def save(n_clicks):
    if n_clicks != None:
        # n = "config.graph = ["
        # for i in [json_serializable(i) for i in config.graph]:
        #     n += i+","
        # n += "]"
        # config.adata.uns["sc_interactive"]["__graph__"] = n

        if "__graph__" not in os.listdir("."):

            os.mkdir("__graph__")

        image_unselected(config.selected)
        file = make_graph_path(config.file_path)
        with open(file,"w") as outfile:
            json_object = json.dumps(config.graph, indent=4, cls=NpEncoder)
            outfile.write(json_object)
        image_selected(config.selected)

        file = make_file_path(config.file_path)
        config.adata.write(file)

    return ""