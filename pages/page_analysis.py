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
from dash.exceptions import PreventUpdate

from app import app

from .graph import *
from .arguments import *

for i in os.listdir("./pages/methods"):
    if i.endswith(".py"):
        f = f"from .methods.{i[:-3]} import *"
        exec(f)

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
                            value=[i for i in config.methods.keys()][0],
                            options=[i for i in config.methods.keys()],
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

def graph2table():
    return [{"Name":i['data']['id'],"Type":i['data']['type'],"Method":i['data']['method']} for i in get_nodes() if i['data']['id'] != 'Raw']

def args2summary(data):

    summary = ""
    for i in data:
        if 'summary' in i.keys():
            summary += f"{i['name']}:{str(i['value'])}"

    return summary

#Add new parameter
methods_implemented = []
for method in config.methods.keys():

    for i in deepcopy(config.methods[method]["args"]["execution"])+\
                deepcopy(config.methods[method]["args"]["postexecution"])+\
                deepcopy(config.methods[method]["args"]["plot"]):

        m_i = (i['name'],i['input'])

        if m_i not in methods_implemented:

            methods_implemented.append(m_i)

            if i['input'] == 'Input':
                input = f"dash.Input('analysis_{i['name']}','value'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = ""
                up_plot = ""
                args = "data"
            elif i['input'] == 'Dropdown':
                input = f"dash.Input('analysis_{i['name']}','value'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = ""
                up_plot = ""
                args = "data"
            elif i['input'] == 'BooleanSwitch':
                input = f"dash.Input('analysis_{i['name']}','on'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = ""
                up_plot = ""
                args = "data"
            elif i['input'] == 'AgTable':
                input = f"dash.Input('analysis_{i['name']}','cellValueChanged'), dash.State('analysis_{i['name']}','rowData'),"
                up = f"if cell != None: data[cell['rowIndex']] = cell['data']; config.active_node_parameters['{i['name']}'] = data"
                # up = "print(cell,data)"
                up_post = "if cell != None: data[cell['rowIndex']] = cell['data']"
                up_plot = up_post
                args = "cell, data"

#Execution
            add_function = f"""
@app.callback(
    dash.Output("dumb","children", allow_duplicate=True),
    {input}
    prevent_initial_call=True
)
def change_parameter_{i['name']}({args}):

    method = get_node(config.selected)["data"]["method"]
    if '{i['name']}' not in [i['name'] for i in config.methods[method]["args"]["execution"]]:
        raise PreventUpdate()

    {up}
    # recomputeWhenUpdate({i['name']},"execution")

    return ""
"""

            exec(add_function, globals(), locals())

#Postexecution and Plot
            add_function = f"""
@app.callback(
    dash.Output("analysis_postargs","children", allow_duplicate=True),
    dash.Output("analysis_plotargs","children", allow_duplicate=True),
    dash.Output("analysis_plot","children", allow_duplicate=True),
    {input}
    prevent_initial_call=True
)
def change_postparameter_{i['name']}({args}):

    if '{i['name']}' in config.block_callback.keys():
        if config.block_callback['{i['name']}']:
            config.block_callback['{i['name']}'] = False
            raise PreventUpdate()

    method = get_node(config.selected)["data"]["method"]
    if '{i['name']}' in [i['name'] for i in config.methods[method]["args"]["postexecution"]]:

        {up_post}
        set_parameters_adata(dict({i['name']}=data))
        set_parameters_node(dict({i['name']}=data))
    
    elif '{i['name']}' in [i['name'] for i in config.methods[method]["args"]["plot"]]:
    
        {up_plot}
        set_plot_parameters_adata(dict({i['name']}=data))
        set_plot_parameters_node(dict({i['name']}=data))

    else:
    
        raise PreventUpdate()

    method = get_node(config.selected)["data"]["method"]

    post_args = deepcopy(config.methods[method]["args"]["postexecution"])
    deactivate = False
    for i in post_args:
        config.block_callback[i['name']] = True
    config.block_callback['{i['name']}'] = False
    # recomputeWhenUpdate({i['name']},"postexecution")
    post_args_object = make_arguments(method, post_args, loaded_args=config.adata.uns[config.selected]['parameters'], add_execution_button=False, add_header="postargs")

    plot_args = deepcopy(config.methods[method]["args"]["plot"])
    deactivate = False
    for i in plot_args:
        config.block_callback[i['name']] = True
    config.block_callback['{i['name']}'] = False
    # recomputeWhenUpdate({i['name']},"plot")
    plot_args_object = make_arguments(method, plot_args, loaded_args=config.adata.uns[config.selected]['plot'], add_execution_button=False, add_header="plot")

    plot = config.methods[method]['plot']()

    return post_args_object, plot_args_object, plot
"""

            exec(add_function, globals(), locals())

#Delete row ag
            if i['input'] == "AgTable" and 'deleteRows' in i.keys():

                if i['deleteRows']:

                    add_function = f"""
@app.callback(
    dash.Output('analysis_{i['name']}', "rowData", allow_duplicate=True),
    dash.Input('analysis_{i['name']}', "virtualRowData"),
    dash.State('analysis_{i['name']}', "rowData"),
    prevent_initial_call=True,
)
def delete_rows_{i['name']}(row, row2):

    if row != None:
        
        config.active_node_parameters["{i['name']}"] = row

    return row
"""
            
                    exec(add_function, globals(), locals())

#Add row ag
            if i['input'] == "AgTable" and 'addRows' in i.keys():

                add_function = f"""
@app.callback(
    dash.Output("analysis_{i['name']}","rowData", allow_duplicate=True),
    dash.Input("analysis_{i['name']}_button","n_clicks"),
    dash.State("analysis_{i['name']}","rowData"),
    prevent_initial_call = True
)
def add_row_{i['name']}(n_clicks, c):

    if n_clicks:

        c.append(
            {i['addRows']}
        )        

    return c
"""
                exec(add_function, globals(), locals())

def set_parameters_adata(args):

    if config.selected not in config.adata.uns.keys():
        config.adata.uns[config.selected] = {"parameters":{}}

    if "parameters" not in config.adata.uns[config.selected].keys():
        config.adata.uns[config.selected]["parameters"] = args
    else:
        for i,j in args.items():
            config.adata.uns[config.selected]["parameters"][i] = j

    config.adata.uns[config.selected]["scinteractive"] = True
    config.adata.uns[config.selected]["method"] = get_node(config.selected)["data"]["method"]

def set_plot_parameters_adata(args):

    if config.selected not in config.adata.uns.keys():
        config.adata.uns[config.selected] = {"plot":{}}

    if "plot" not in config.adata.uns[config.selected].keys():
        config.adata.uns[config.selected]["plot"] = args
    else:
        for i,j in args.items():
            config.adata.uns[config.selected]["plot"][i] = j

    config.adata.uns[config.selected]["scinteractive"] = True
    config.adata.uns[config.selected]["method"] = get_node(config.selected)["data"]["method"]

def set_parameters_node(args):

    pos = get_node_pos(config.selected)

    for i,j in args.items():
        config.graph[pos]["data"]["parameters"][i] = j

def set_plot_parameters_node(args):

    pos = get_node_pos(config.selected)

    for i,j in args.items():
        config.graph[pos]["data"]["plot"][i] = j

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

    args = deepcopy(config.methods[method]["args"]["execution"])
    args = args_eval(args)

    h5ad_file = get_node(config.selected)["data"]["h5ad_file"]

    node = {
        'data': {
            'id': name, 
            'name': method,
            'method':method, 
            'type':config.methods[method]["properties"]['type'], 
            'color': "blue",
            'h5ad_file':h5ad_file,
            'image': '',
            'computed':False,
            'opacity':.3,
            'summary':'', 
            'parameters':args,
            'plot':{}
        }, 
        'position':{'x':config.max_x + 30,'y':0},
    }
    #make active node that is the one that will be presented
    # nodes_update_pos(config.graph,cytoscape)

    config.active_node_parameters = args.copy()
    config.graph.append(node)

    l = load_node(name)
    # make_node_summary(config.selected)

    return config.graph, graph2table(), l

#Execute analysis button
@app.callback(
    dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
    dash.Output('execute-input-modal', 'is_open', allow_duplicate=True),
    dash.Output('execute-computed-modal', 'is_open', allow_duplicate=True),
    dash.Output('analysis_postargs', 'children', allow_duplicate=True),
    dash.Output('analysis_plotargs', 'children', allow_duplicate=True),
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

                #Get incoming input, old input and method
                input = config.active_node_parameters["input"]
                old_input = get_node(config.selected)["data"]["parameters"]["input"]
                method = get_node(config.selected)["data"]["method"]

                #Get erguments of incoming
                inputArgs = get_args(input)

                #Execute code
                outputArgs = config.methods[method]["function"](config.adata, inputArgs, config.active_node_parameters)

                # Remove old edge if any and add new
                edge_rm(old_input, config.selected)
                edge_add(input, config.selected)

                #save parameters
                set_output(outputArgs)
                set_parameters_adata(config.active_node_parameters)
                set_parameters_node(config.active_node_parameters)

                #Add post arguments
                post_args = deepcopy(config.methods[method]["args"]["postexecution"])
                post_args = args_eval(post_args)
                set_parameters_adata(post_args)
                set_parameters_node(post_args)

                #Add plot arguments
                plot_args = deepcopy(config.methods[method]["args"]["plot"])
                plot_args = args_eval(plot_args)
                set_plot_parameters_adata(plot_args)
                set_plot_parameters_node(plot_args)

                #Activate node
                activate_node(config.selected)
                deactivate_downstream(config.selected)

        post_args_object = {}
        plot = []
        node_data = get_node(config.selected)['data']
        if node_data['computed']:
            post_args = deepcopy(config.methods[node_data['method']]["args"]["postexecution"])
            plot_args = deepcopy(config.methods[node_data['method']]["args"]["plot"])
            config.block_callback = {}
            for i in post_args:
                config.block_callback[i["name"]] = True
            for i in plot_args:
                config.block_callback[i["name"]] = True

            post_args_object = make_arguments(node_data['method'], post_args, add_execution_button=False, add_header="postargs")
            plot_args_object = make_arguments(node_data['method'], plot_args, add_execution_button=False, add_header="plot")

            plot = config.methods[node_data['method']]["plot"]()

    else:

        raise PreventUpdate()

    return config.graph, warning_input, warning_computed, post_args_object, plot_args_object, plot

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

#Load button
@app.callback(
    dash.Output('graph_cytoscape', 'elements', allow_duplicate=True),
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

        node_data = get_node(config.selected)['data']
        if node_data['computed']:
            post_args = deepcopy(config.methods[node_data['method']]["args"]["postexecution"])
            plot_args = deepcopy(config.methods[node_data['method']]["args"]["plot"])
            config.block_callback = {}
            for i in post_args:
                config.block_callback[i["name"]] = False
            for i in plot_args:
                config.block_callback[i["name"]] = False

    return config.graph, l

@app.callback(
    dash.Output('graph_dropdown_load', 'value', allow_duplicate=True),
    dash.Input('graph_cytoscape', 'tapNode'),
    prevent_initial_call=True
)
def display_click_data(tap_node_data):

    if tap_node_data is not None:
        return tap_node_data['data']['id']
    
    return None

#Save analysis
@app.callback(
    dash.Output("dumb","children",allow_duplicate=True),
    dash.Input("graph_save_button","n_clicks"),
    prevent_initial_call=True
)
def save(n_clicks):
    if n_clicks != None:

        save_graph()

        adapt_adata_saving()
        save_adata()
        adapt_adata_loaded()

    return ""

@app.callback(
    dash.Output("dumb","children",allow_duplicate=True),
    dash.Input("analysis_saveimage_button","n_clicks"),
    dash.State("analysis_plot","children"),
    prevent_initial_call=True
)
def savefigure(n_clicks,fig):

    if n_clicks != None:

        figs = get_figures(fig)

        count = 0
        for fig in figs:
            namefig = f"{config.analysis_folder}/report/figures/{config.selected}_{count}.png"
            while os.path.isfile(namefig):
                count += 1
                namefig = f"{config.analysis_folder}/report/figures/{config.selected}_{count}.png"
            fig.write_image(namefig)

            config.report += f"![](./figures/{config.selected}_{count}.png)\n"

    return ""
