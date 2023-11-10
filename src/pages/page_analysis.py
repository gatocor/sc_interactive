from dash import dcc, html, dash_table, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_editor_components

import os
from traceback import format_exc
import base64
from io import BytesIO
from matplotlib.figure import Figure

from app import app
from app import *

from general import *

methods_menus = []
d = {}
for folder, dirs, files in os.walk("./methods"):

    if not folder.endswith("__pycache__"):
        d[folder] = []

        for i in files:
            if i.endswith(".py") and i != "__init__.py":
                module = f"{folder.replace('.py','').replace('./','').replace('/','.')}.{i.replace('.py','')}"
                f = f"from {module} import *"
                exec(f)

                id = folder[1:]+"/"+i.replace('.py','')
                d[folder].append(dbc.DropdownMenuItem(i.replace('.py',''),id=id))
                methods_menus.append(id)

        folder_parent = folder
        while len(folder_parent.split("/")) > 2:
            name = ""
            for i in folder_parent.split("/")[:-1]:
                name += i+"/"
            folder_kid = folder_parent.split("/")[-1]
            folder_parent = name[:-1]

            d[folder_parent].append(dbc.DropdownMenu(d[folder],label=folder_kid,size="md",className="mb-3",color="light"))

for folder, dirs, files in os.walk("./methods_plot"):

    for i in files:
        if i.endswith(".py") and i != "__init__.py":
            module = f"{folder.replace('.py','').replace('./','').replace('/','.')}.{i.replace('.py','')}"
            f = f"from {module} import *"
            exec(f)

############################################################################################################################################
############################################################################################################################################
# Layout
############################################################################################################################################
############################################################################################################################################
tab_algorithm = dbc.Card(
    [
        dbc.CardBody(
            id='analysis_args',
            children = [],
            # width='50%',
        ),
        dbc.CardFooter(
            children=[dbc.Button("Execute",id="analysis_execute_button")],
        ),
    ],
    color="#CED4DA",
    className="mt-3",
)
tab_plot = dbc.Row([
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    dcc.Dropdown(
                        id='analysis_plot_dropdown',
                        value = None,
                        clearable = True,
                        options = []
                    )
                    # width='50%',
                ),
                dbc.CardBody(
                    id='analysis_plotargs',
                    children = [],
                    # width='50%',
                ),
                dbc.CardFooter(
                    children=[
                        dbc.Button("Execute",id="analysis_execute_plot_button"),
                    ],
                ),
            ],
            color="#CED4DA",
            className="mt-3",
        ),
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    id='analysis_plot',
                    children = [],
                    # width='50%',
                ),
                dbc.CardFooter(
                    children=[
                        dbc.Button("Save Image",id="analysis_saveimage_button"),
                    ],
                ),
            ],
            color="#CED4DA",
            className="mt-3",
        )
    )
])
tab_saved_images = dbc.Card(
    [
        dbc.CardBody(
            id='analysis_plots_list',
            children = [],
            # width='50%',
        ),
        dbc.CardFooter(
            children=[
                dcc.Dropdown(id="analysis_saved_plots_dropdown",options=[]),
                dbc.Button("Delete Image",id="analysis_saved_plots_remove",style={"background-color":"red"}),
                dbc.Button("Reload Image",id="analysis_saved_plots_upload"),
            ],
        ),
    ],
    color="#CED4DA",
    className="mt-3",
)
tab_report = dbc.Card(
    children = [
        dbc.CardHeader(
            children = dbc.Switch(id="analysis-show-editor",value=False,label="Activate editor"),
        ),
        dbc.CardBody(
            id='analysis_report',
            children = [],
            # width='50%',
        ),
        dbc.CardFooter(
            children = dbc.Button(id="analysis-save-report",children="Save Report",disabled=False),
        ),
    ],
    color="#CED4DA",
    className="mt-3",
)
tab_h5ad = dbc.Card(
    [
        dbc.CardBody(
            id='analysis_inspector',
            children = [],
            # width='50%',
        )
    ],
    color="#CED4DA",
    className="mt-3",
)
tab_info = dbc.Card(
    [
        dbc.CardBody(
            id='analysis_info',
            children = [],
            # width='50%',
        )
    ],
    color="#CED4DA",
    className="mt-3",
)
tab_info_plot = dbc.Card(
    [
        dbc.CardBody(
            id='analysis_info_plot',
            children = [],
            # width='50%',
        )
    ],
    color="#CED4DA",
    className="mt-3",
)

def layout():
    return dbc.Container(
        [
            html.P(id="dumb",children=""),
            dbc.Row(
                children=[
                    dbc.Col(html.H1(id="title",children="Graphical Analysis")),
                    dbc.Col(
                        dbc.Button(id='graph_save_button', children="Save",
                                    size="lg",
                                    style={
                                        "background-color":"#343A40",
                                        'width': '280px', 
                                    }      
                                    ),
                        width={"offset":5}
                    )
                ],
                justify="center",
                className="mb-4"
            ),
            dbc.Modal(
                [
                    dbc.ModalHeader(close_button=False,
                        children=[html.Div("New"), dbc.Switch(id="new_h5ad",value=False,label="New h5ad")]             
                    ),
                    dbc.ModalBody(id="new-message",
                                children=[
                                    html.Div("You are creating a new node. Please select a node to which attach this one and press add."),
                                    html.Label("Source node"),
                                    dcc.Dropdown(id="new-dropdown-before",value="Raw",options=node_names(), clearable=False),
                                    html.Label("Target node (can be None for a new branch)"),
                                    dcc.Dropdown(id="new-dropdown-after",value="Raw",options=node_names(), clearable=True)
                                ]
                    ),
                    dbc.ModalFooter([
                            dbc.Button("New", id="new-proceed", className="ml-auto"),
                            dbc.Button("Cancel", id="new-cancel", className="ml-auto"),
                    ])
                ],
                backdrop=False,
                id="new-modal",
                size="sm",
            ),
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
                    dbc.ModalBody(children=html.Div("The previous cells must be computed before computing this one.")),
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
                children=dbc.ButtonGroup([
                        dcc.Dropdown(
                            id='graph_dropdown_load',
                            value=None,
                            options=node_names(),
                            clearable=True,
                            className="mb-3"
                        ),
                        dbc.Button(
                            id='graph_load_button',
                            children="Load",
                            className="mb-3"
                        ),
                        dbc.Button(
                            children="Delete",
                            id="analysis_delete_button", 
                            style={"background-color":"red"}, 
                            className="mb-3"
                        ),
                        dbc.DropdownMenu([d["./methods"][0]],label="Add New Method",id = 'graph_dropdown_analysis',className="mb-3",),
                ],size="md",className="mb-2")
            ),
            dbc.Row(
                cyto.Cytoscape(
                    id='graph_cytoscape',
                    # layout={'name':'breadthfirst','roots':'[id="Raw"]'},#{'name': 'preset'},
                    layout=GRAPHLAYOUT,#{'name': 'preset'},
                    style={'width': '100%', 'height': '700px'},
                    elements=config.graph,
                    userZoomingEnabled=False,  # Disable zooming
                    userPanningEnabled=True,  # Disable panning
                    stylesheet=GRAPHSTYLESHEET,
                ),
                style={'width': '100wh', 'height':'120wh', 'margin': 'auto'}, #vh/wh
            ),
            dbc.Alert(
                children="",
                id="execution-error",
                dismissable=True,
                is_open=False,
                color="danger"
            ),
            dbc.Alert(
                id="node_analysis",
                children=[
                    html.H1(id="analysis_name",children=["Raw"]),
                    dbc.Tabs(
                        [
                            dbc.Tab(tab_algorithm, label="Algorithm"),
                            dbc.Tab(tab_plot, label="Plot"),
                            dbc.Tab(tab_saved_images, label="Saved Plots"),
                            dbc.Tab(tab_report, label="Report"),
                            dbc.Tab(tab_h5ad, label="Object Inspector"),
                            dbc.Tab(tab_info, label="Algorithm Information"),
                            dbc.Tab(tab_info_plot, label="Plot Algorithm Information"),
                        ]
                    ),
                ],
                color="dark", 
                is_open=False
            ),
        ],
        fluid=True
    )

############################################################################################################################################
############################################################################################################################################
# Functions
############################################################################################################################################
############################################################################################################################################
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

def parameters_eval(args, populate = None, count = 0):
    
    if count == 0:
        args = deepcopy(args)

    if type(args) == list:
        for i,j in enumerate(args):
            args[i] = parameters_eval(j, populate, count+1)

            if count == 0:
                value = get_value(args[i])

                if populate == "execution":
                    config.active_node_parameters[args[i]["name"]] = value
                elif populate == "plot":
                    config.active_plot_parameters[args[i]["name"]] = value

    elif type(args) == dict:
        if "function" in args.keys():
            return fvalue(args)
        else:
            for i,j in args.items():
                args[i] = parameters_eval(j, populate, count+1)

    return args

def parameters2args(args, populate):

    args = parameters_eval(args, populate)

    return {i["name"]:get_value(i) for i in args}    

def method_create_pars(args_type):

    if args_type == "execution":
        method_args = config.methods[get_node(config.selected)["data"]["method"]]["args"]
        for j in method_args:
            config.active_node_parameters[j["name"]] = j["properties"]["value"]
    elif args_type == "plot":
        if config.selected_plot != None:
            method_args = config.methods_plot[config.selected_plot]["args"]
            for j in method_args:
                config.active_plot_parameters[j["name"]] = j["properties"]["value"]
            else:
                return None

    args = parameters_eval(method_args, args_type)

def make_arguments(id, arg_list, loaded_args={}, plot=False, add_header="args"):

    if arg_list == []:
        return []

    #Avoid updates after creating the inputs
    arglist = parameters_eval(arg_list)
    for i in arglist:
        config.block_callback[i["name"]] = True

    l = []

    for i,arg in enumerate(arglist):

        if fvisible(arg):
        
            l.append(
                dbc.Tooltip(
                        arg["description"],
                        target=id+str(i),
                        placement="bottom",
                )
            )
            lab = html.Label(arg["name"],id=id+str(i),style={'text-align': 'right'})
            if arg["input"] == "Input":

                try:
                    arg["properties"]["value"] = loaded_args[arg["name"]]
                except:
                    loaded_args[arg["name"]] = arg["properties"]["value"]

                input = dbc.Input(
                            id="analysis_"+str(arg["name"]),
                            **arg["properties"]
                        )

            # elif arg["input"] == "Dropdown":

            #     try:
            #         arg["properties"]["value"] = loaded_args[arg["name"]]
            #     except:
            #         loaded_args[arg["name"]] = arg["properties"]["value"]

            #     input = dcc.Dropdown(
            #                 id="analysis_"+str(arg["name"]),
            #                 **arg["properties"]
            #             )
                
            # elif arg["input"] == "BooleanSwitch":

            #     try:
            #         arg["properties"]["on"] = loaded_args[arg["name"]]
            #     except:
            #         loaded_args[arg["name"]] = arg["properties"]["on"]

            #     input = daq.BooleanSwitch(
            #                 id="analysis_"+str(arg["name"]),
            #                 **arg["properties"]
            #             )

            # elif arg["input"] == "AgTable":

            #     try:
            #         arg["properties"]["data"] = loaded_args[arg["name"]]
            #     except:
            #         None

            #     if "deleteRows" in arg.keys():
            #         d = arg["deleteRows"]
            #     else:
            #         d = False

            #     input = [
            #         dbc.Row(
            #             html.Div(
            #                     children= ag_table("analysis_"+str(arg["name"]), arg["properties"]["header"], arg["properties"]["data"], deleteRows=d),
            #             ),
            #         ),
            #     ]

            #     if "addRows" in arg.keys():

            #         input += [
            #             dbc.Row(
            #                 html.Div(
            #                         children= dbc.Button(id="analysis_"+str(arg["name"])+"_button", children="Add"),
            #                 ),
            #             ),
            #         ]

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

    if not plot:
        if config.show_parameters:
            l.append(
                html.Div(
                    dbc.Button("Show Less Arguments",id="analysis_unfold_execution_button",color="dark"),
                    className="d-grid gap-2"
                )
            )
        else:
            l.append(
                html.Div(
                    dbc.Button("Show More Arguments",id="analysis_unfold_execution_button",color="dark"),
                    className="d-grid gap-2"
                )
            )
    else:
        if config.show_plot:
            l.append(
                html.Div(
                    dbc.Button("Show Less Arguments",id="analysis_unfold_plot_button",color="dark"),
                    className="d-grid gap-2"
                )
            )
        else:
            l.append(
                html.Div(
                    dbc.Button("Show More Arguments",id="analysis_unfold_plot_button",color="dark"),
                    className="d-grid gap-2"
                )
            )

    return l

def load_node(name):

    if name == "Raw":
        return [], [], []

    method = get_node(name)['data']['method']
    args = get_node(name)['data']['parameters']
    l = make_arguments(method, config.methods[method]["args"], args)
    l3 = []
    p = []
            
    return l, l3, p

############################################################################################################################################
############################################################################################################################################
# Callbacks
############################################################################################################################################
############################################################################################################################################

#Add new parameter
methods_implemented = []
for method in config.methods.keys():

    for i in deepcopy(config.methods[method]["args"]):

        m_i = (i['name'],i['input'])

        if m_i not in methods_implemented:

            methods_implemented.append(m_i)

            if i['input'] == 'Input':
                input = f"Input('analysis_{i['name']}','value'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                args = "data"

            add_input = ""
            if i["name"] == "input":
                add_input = f"config.h5ad_file = None"

#Execution
            add_function = f"""
@app.callback(
    Output("analysis_args","children", allow_duplicate=True),
    {input}
    prevent_initial_call=True
)
def change_parameter_{i['name']}({args}):

    method = get_node(config.selected)["data"]["method"]
    if '{i['name']}' not in [i['name'] for i in config.methods[method]["args"]]:
        raise PreventUpdate()

    if config.block_callback['{i['name']}']:
        config.block_callback['{i['name']}'] = False
        raise PreventUpdate()

    {up}
    clean_arguments('{i['name']}', config.methods[method]["args"], config.active_node_parameters)
    args_object = make_arguments(method, config.methods[method]["args"], config.active_node_parameters)
    config.block_callback['{i['name']}'] = False

    return args_object
"""

            exec(add_function, globals(), locals())

methods_plot_implemented = []
for method in config.methods_plot.keys():

    for i in deepcopy(config.methods_plot[method]["args"]):

        m_i = (i['name'],i['input'])

        if m_i not in methods_plot_implemented:

            methods_plot_implemented.append(m_i)

            if i['input'] == 'Input':
                input = f"Input('analysis_{i['name']}','value'),"
                up = f"config.active_plot_parameters['{i['name']}'] = data"
                args = "data"

# Plot
            add_function = f"""
@app.callback(
    # Output("analysis_plotargs","children", allow_duplicate=True),
    Output("dumb","children", allow_duplicate=True),
    {input}
    prevent_initial_call=True
)
def change_parameter_{i['name']}({args}):

    # method = get_node(config.selected)["data"]["method"]
    # if '{i['name']}' not in [i['name'] for i in config.methods[method]["args"]]:
    #     raise PreventUpdate()

    # if config.block_callback['{i['name']}']:
    #     config.block_callback['{i['name']}'] = False
    #     raise PreventUpdate()

    {up}
    # clean_arguments('{i['name']}', config.methods_plot[method]["args"], config.active_plot_parameters)
    # args_object = make_arguments(method, config.methods_plot[method]["args"], config.active_plot_parameters)
    # config.block_callback['{i['name']}'] = False

    return ""#args_object
"""

            exec(add_function, globals(), locals())

# #Delete row ag
#             if i['input'] == "AgTable" and 'deleteRows' in i.keys():

#                 if i['deleteRows']:

#                     add_function = f"""
# @app.callback(
#     Output('analysis_{i['name']}', "rowData", allow_duplicate=True),
#     Input('analysis_{i['name']}', "virtualRowData"),
#     State('analysis_{i['name']}', "rowData"),
#     prevent_initial_call=True,
# )
# def delete_rows_{i['name']}(row, row2):

#     if row != None:
        
#         config.active_node_parameters["{i['name']}"] = row

#     return row
# """
            
#                     exec(add_function, globals(), locals())

# #Add row ag
#             if i['input'] == "AgTable" and 'addRows' in i.keys():

#                 add_function = f"""
# @app.callback(
#     Output("analysis_{i['name']}","rowData", allow_duplicate=True),
#     Input("analysis_{i['name']}_button","n_clicks"),
#     State("analysis_{i['name']}","rowData"),
#     prevent_initial_call = True
# )
# def add_row_{i['name']}(n_clicks, c):

#     if n_clicks:

#         c.append(
#             {i['addRows']}
#         )        

#     return c
# """
#                 exec(add_function, globals(), locals())

#Update Dropdown (change of cytoscape)
@app.callback(
    Output('graph_dropdown_load', 'options', allow_duplicate=True),
    [
        Input('graph_cytoscape', 'elements')
    ],
    prevent_initial_call=True
)
def graph_update_dropdown(_):

    return node_names()

#Update Table (change of cytoscape)
@app.callback(
    Output('graph_table', 'data', allow_duplicate=True),
    [
        Input('graph_cytoscape', 'elements')
    ],
    prevent_initial_call=True
)
def graph_update_table(_):

    return graph2table()

@app.callback(
    Output('dumb', 'children', allow_duplicate=True),
    Input('graph_cytoscape', 'tapNode'),
    prevent_initial_call=True
)
def graph_update_pos(element):

    node_update_pos(config.graph, element['data']['id'], element['position'])
    config.max_x = max([i['position']['x'] for i in get_nodes()])

    return ""

#Add node
@app.callback(
    Output('new-dropdown-after', 'value', allow_duplicate=True),
    Output('new-dropdown-after', 'options', allow_duplicate=True),
    Input('new-dropdown-before', 'value'),
    prevent_initial_call=True
)
def graph_new_modal(innode):

    l = [i["data"]["id"] for i in get_node_children(innode)]

    return None, l

code = ""
for i in methods_menus:

    code += f"""
    Input("{i}", "n_clicks"),"""

code = f"""@app.callback(
    Output('new-modal', 'is_open', allow_duplicate=True),
    Output('new-dropdown-before', 'value', allow_duplicate=True),
    Output('new-dropdown-before', 'options', allow_duplicate=True),
    # Output("graph_dropdown_analysis", "children"),
    {code}
    State("graph_dropdown_load","value"),
    prevent_initial_call=True
)
def graph_new_node(*_):
    ctx = dash.callback_context
    load = ctx.triggered[0]["prop_id"].split(".")[0].split("/")[-1]

    print(ctx.triggered[0]["prop_id"])
    config.add_method = load

    return True, _[-1], node_names()
"""
# print(code)
exec(code)

# code = ""
# for i in methods_menus:

#     code += f"""
#     Input("{i}", "n_clicks"),"""

# code = f"""@app.callback(
#     Output("graph_dropdown_analysis", "isOpen"),
#     # Output("graph_dropdown_analysis", "children"),
#     {code}
#     prevent_initial_call=True
# )
# def hide_show_submenu(*_):
#     ctx = dash.callback_context
#     input_id = ctx.triggered[0]["prop_id"].split(".")[0]
#     # return "d-block" if input_id == "graph_dropdown_analysis" else "d-none"
#     return False
# """
# print(code)
# exec(code)

@app.callback(
    Output('new-modal', 'is_open', allow_duplicate=True),
    Input('new-cancel', 'n_clicks'),
    prevent_initial_call=True
)
def graph_new_node(_):

    return False

@app.callback(
    Output('node_analysis', 'is_open', allow_duplicate=True),
    Output('new-modal', 'is_open', allow_duplicate=True),
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('graph_table', 'data', allow_duplicate=True),
    Output('analysis_name','children',allow_duplicate=True),
    Output('analysis_args', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    Output('analysis_inspector', 'children', allow_duplicate=True),
    Output('analysis_plot_dropdown', 'options', allow_duplicate=True),
    Output('analysis_plot_dropdown', 'value', allow_duplicate=True),
    Output('new_h5ad', 'value', allow_duplicate=True),
    [
        Input('new-proceed', 'n_clicks')
    ],
    [
        State('new-dropdown-before', 'value'),
        State('new-dropdown-after', 'value'),
        State('new_h5ad', 'value'),
        State('graph_cytoscape', 'layout'),
    ],
    prevent_initial_call=True
)
def graph_new_node(_, input, output, new, cytoscape):

    #Make new name
    method = config.add_method
    name = config.add_method
    nodes = node_names()
    count = 0
    while name in nodes:
        count += 1
        name = f"{method}_{str(count)}"

    if new:
        typ = "QC_saved"
    else:
        typ = "QC"

    node = {
        'data': {
            'id': name, 
            'name': method,
            'method':method, 
            'type':typ, 
            'color': "blue",
            'h5ad_file':None,
            'image': '',
            'computed':False,
            'opacity':.3,
            'summary':'', 
            'parameters':{"input":input},
            'plots':[],
            'report':"",
            'new_h5ad':new
        }, 
        'position':{'x':config.max_x + 30,'y':0},
    }
    #make active node that is the one that will be presented
    # nodes_update_pos(config.graph,cytoscape)

    if input != config.h5ad_file: #Load appropiate node
        # save_adata()
        innode = get_node(input)["data"]
        load_adata(f"{config.analysis_folder}/h5ad/{innode['h5ad_file']}")

    config.report = ""
    #Create node
    config.graph.append(node)
    #Create block callbacks
    config.block_callback = {}
    #Selected
    unselect_node(config.selected)
    config.selected = name
    select_node(config.selected)
    #Create parameters
    method_create_pars("execution")
    #Create active
    config.active_node_parameters["input"] = input
    set_parameters(config.active_node_parameters, "parameters")
    # set_active_node_parameters({i["name"]:get_value(i) for i in args.copy()})
    #Assign new file
    innode = get_node(input)
    pos = get_node_pos(config.selected)
    if node["data"]["new_h5ad"]:
        config.graph[pos]["data"]["h5ad_file"] = new_h5ad_file()
    elif innode["data"]["id"] == "Raw":
        config.graph[pos]["data"]["h5ad_file"] = new_h5ad_file()
    else:
        config.graph[pos]["data"]["h5ad_file"] = innode["data"]["h5ad_file"]
    
    edge_add(input, config.selected)

    if output:

        deactivate_downstream(output)
        deactivate_node(output)
        node_reassign_input(output,name)

    #load interactive plots
    l, l3, p = load_node(name)

    inspector = print_to_string(config.adata)

    plot_options = get_plot_methods(method)

    return True, False, config.graph, graph2table(), name, l, l3, p, html.Pre(inspector, style={"white-space":"pre-wrap"}), plot_options, None, False

#Execute analysis button
@app.callback(
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('execute-computed-modal', 'is_open', allow_duplicate=True),
    Output('execution-error', 'is_open', allow_duplicate=True),
    Output('execution-error', 'children', allow_duplicate=True),
    [
        Input('analysis_execute_button', 'n_clicks')
    ],
    State('execute-computed-modal', 'is_open'),
    prevent_initial_call=True
)
def execute(n_clicks, warning_computed):

    if n_clicks != None:

        for node in get_node_ancestors(config.selected):
            
            if not node['data']['computed']:
            
                return config.graph, True, False, []
            
        #Get incoming input, old input and method
        input = config.active_node_parameters["input"]

        method = get_node(config.selected)["data"]["method"]

        #Execute code
        innode = get_node(input)["data"]
        pos = get_node_pos(config.selected)

        if innode["h5ad_file"] != config.h5ad_file: #Load appropiate node
            load_adata(f"{config.analysis_folder}/h5ad/{innode['h5ad_file']}")
        config.h5ad_file = config.graph[pos]["data"]["h5ad_file"]      

        try:
            config.methods[method]["function"](config.adata, config.active_node_parameters)
        except BaseException as e:
            l = [html.H2("Node failed to execute.")]
            for i in format_exc().split("\n"):
                l.append(html.Pre(i))
            return config.graph, warning_computed, True, l

        #save parameters
        set_parameters(config.active_node_parameters, 'parameters')

        #Activate node
        activate_node(config.selected)
        deactivate_downstream(config.selected)

        save_adata()
        save_graph()

        return config.graph, warning_computed, False, ""

    else:

        raise PreventUpdate()

@app.callback(
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Input('analysis_plot_dropdown', 'value'),
    prevent_initial_call=True
)
def graph_new_node(plot_type):

    if plot_type != None:

        #Create parameters
        config.selected_plot = plot_type

        method_create_pars("plot")

        l = make_arguments(method, config.methods_plot[plot_type]["args"], config.active_plot_parameters, plot=True)

        return l
    else:
        raise PreventUpdate()

@app.callback(
    Output('analysis_plot', 'children', allow_duplicate=True),
    Output('execution-error', 'is_open', allow_duplicate=True),
    Output('execution-error', 'children', allow_duplicate=True),
    Input('analysis_execute_plot_button', 'n_clicks'),
    State('analysis_plot_dropdown','value'),
    prevent_initial_call=True
)
def execute_plot(n_clicks, plot_type):

    if n_clicks != None:
            
        try:
            fig = config.methods_plot[plot_type]['function']()
            config.figure = fig

            if isinstance(fig,Figure):
                # Save it to a temporary buffer.
                buf = BytesIO()
                fig.savefig(buf, format="png", transparent=True)
                # Embed the result in the html output.
                fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
                fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
                plot =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)
            else:
                print("bonico")

        except BaseException as e:
            l = [html.H2("Node failed to execute.")]
            for i in format_exc().split("\n"):
                l.append(html.Pre(i))
            return [], True, l

        return plot, False, []

    else:

        raise PreventUpdate()

#Delete analysis button
@app.callback(
    Output('delete-modal', 'is_open'),
    Output('delete-message', 'children'),
    Input('analysis_delete_button', 'n_clicks'),
    State('graph_dropdown_load', 'value'),
    prevent_initial_call=True
)
def delete(n_clicks, value):
        
    if n_clicks != None:

        l = [
                html.Div("You are about to remove an analysis node. This will remove all the computations performed over this node and descending nodes of the analysis. Are you sure that you want to proceed?"),
                html.Div(),
                html.Div(f"Analysis to be removed: {value}")
            ]

        return True, l
    
    return False, []

#Proceed modal delete
@app.callback(
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('analysis_args', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    Output('delete-modal', 'is_open', allow_duplicate=True),
    Output('graph_dropdown_load', 'value', allow_duplicate=True),
    Input('delete-proceed', 'n_clicks'),
    State('graph_dropdown_load', 'value'),
    prevent_initial_call=True
)
def delete_confirmation(n_clicks, val):
    
    if n_clicks != None:

        if val not in  ["Raw", None]:

            deactivate_downstream(val)
            node_rm(val)

            if val == config.selected:
                config.selected = 'Raw'
            
            modal = False

            l, l3, p = load_node(config.selected)

            return config.graph, l, l3, p, modal, None
        
        else:

            l, l3, p = load_node(config.selected)
            modal = False

            return config.graph, l, l3, p, modal, val
    
    else:

        raise PreventUpdate()

#Cancel modal delete
@app.callback(
    Output('delete-modal', 'is_open', allow_duplicate=True),
    [
     Input('delete-cancel', 'n_clicks'),
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
    Output('node_analysis', 'is_open', allow_duplicate=True),
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('analysis_name','children',allow_duplicate=True),
    Output('analysis_args', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    Output('analysis_inspector', 'children', allow_duplicate=True),
    Output('analysis_info', 'children', allow_duplicate=True),
    Output('analysis_plot_dropdown', 'options', allow_duplicate=True),
    Output('analysis_plot_dropdown', 'value', allow_duplicate=True),
    Output("analysis_report","children", allow_duplicate=True),
    Output("analysis_plots_list", "children", allow_duplicate=True),
    [
        Input('graph_load_button', 'n_clicks')
    ],
    [
        State('graph_dropdown_load', 'value'),
    ],
    prevent_initial_call=True
)
def load_analysis(_, name):

    if _ != None and name != None:

        #adata
        node = get_node(name)["data"]
        method = node["method"]
        if node["h5ad_file"] != config.h5ad_file:
            #reload file previous to the model
            load_adata(f"{config.analysis_folder}/h5ad/{node['h5ad_file']}")
            config.h5ad_file = node["h5ad_file"]
            # adapt_adata_loaded()

        #Create block callbacks
        config.block_callback = {}
        #Create active
        pos = get_node_pos(name)
        config.active_node_parameters = deepcopy(config.graph[pos]["data"]["parameters"])
        config.report = deepcopy(config.graph[pos]["data"]["report"])
        #Selected
        unselect_node(config.selected)
        config.selected = name
        select_node(config.selected)

        # #Create parameters
        # method_create_pars(["execution"])
        # # set_active_node_parameters({i["name"]:get_value(i) for i in args.copy()})

        l, l3, p = load_node(name)

        inspector = print_to_string(config.adata)

        plot_options = get_plot_methods(method)

        info = config.methods[method]["docs"]

        report = [
            dbc.Col(id="analysis-markdown", children=markdown_to_dash(config.report)),
        ]

        cols = [
            {"field":"style"},
            {"field":"fig"},
        ]

        table = dag.AgGrid(
                rowData=config.graph[pos]['data']['plots'],
                columnDefs=cols,
                defaultColDef={"editable": True},
                columnSize="sizeToFit",
            )

        return True, config.graph, name, l, l3, p, html.Pre(inspector, style={"white-space":"pre-wrap"}), html.Pre(info, style={"white-space":"pre-wrap"}), plot_options, None, report, table
    
    else:
        raise PreventUpdate()

@app.callback(
    Output('graph_dropdown_load', 'value', allow_duplicate=True),
    Input('graph_cytoscape', 'tapNode'),
    prevent_initial_call=True
)
def display_click_data(tap_node_data):

    if tap_node_data is not None:
        return tap_node_data['data']['id']
    
    return None

#Save analysis
@app.callback(
    Output("dumb","children",allow_duplicate=True),
    Input("graph_save_button","n_clicks"),
    prevent_initial_call=True
)
def save(n_clicks):
    if n_clicks != None:

        save_adata()
        save_graph()
        clean_h5ad()

    return ""

@app.callback(
    Output("analysis-show-editor", "children", allow_duplicate=True),
    Output("analysis_report","children", allow_duplicate=True),
    Output("analysis_plots_list", "children", allow_duplicate=True),
    Input("analysis_saveimage_button","n_clicks"),
    State("analysis_plot","children"),
    State("analysis_plot_dropdown","value"),
    prevent_initial_call=True
)
def savefigure(n_clicks,fig,style):

    if n_clicks != None and style != None:

        if isinstance(config.figure,Figure):

            count = 0
            namefig = f"{config.analysis_folder}/report/figures/{config.selected}_{count}.png"
            while os.path.isfile(namefig):
                count += 1
                namefig = f"{config.analysis_folder}/report/figures/{config.selected}_{count}.png"
            config.figure.savefig(namefig,transparent=True)            

        else:

            figs = get_figures(fig)

            count = 0
            for fig in figs:
                namefig = f"{config.analysis_folder}/report/figures/{config.selected}_{count}.png"
                while os.path.isfile(namefig):
                    count += 1
                    namefig = f"{config.analysis_folder}/report/figures/{config.selected}_{count}.png"
                fig.write_image(namefig)


        pos = get_node_pos(config.selected)
        config.graph[pos]['data']['plots'].append(
            {
                "style":style,
                "parameters":deepcopy(config.active_plot_parameters),
                "fig":f"{config.selected}_{count}.png"
            }
        )
        config.report += f"![](./figures/{config.selected}_{count}.png)\n"
        config.graph[pos]['data']['report'] = config.report

        save_graph()

        report = dbc.Col(id="analysis-markdown", children=markdown_to_dash(config.report))

        table = dag.AgGrid(
                rowData=config.graph[pos]['data']['plots'],
                columnDefs=[{"field":"style"},{"field":"fig"}],
                defaultColDef={"editable": True},
                columnSize="sizeToFit",
            )

        return False, report, table
    
    else:

        raise PreventUpdate()

@app.callback(
    Output("analysis_args","children",allow_duplicate=True),
    Input("analysis_unfold_execution_button","n_clicks"),
    prevent_initial_call=True
)
def unfold_execution(n_clicks):

    if n_clicks != None:

        config.show_parameters = not config.show_parameters

        method = get_node(config.selected)["data"]["method"]
        args_object = make_arguments(method, config.methods[method]["args"], config.active_node_parameters)

        return args_object

    else:

        raise PreventUpdate()

@app.callback(
    Output("analysis_plotargs","children",allow_duplicate=True),
    Input("analysis_unfold_plot_button","n_clicks"),
    prevent_initial_call=True
)
def unfold_plot(n_clicks):

    if n_clicks != None:

        config.show_plot = not config.show_plot

        method = config.selected_plot
        args_object = make_arguments(method, config.methods_plot[method]["args"], config.active_plot_parameters, plot=True)

        return args_object

    else:

        raise PreventUpdate()

@app.callback(
    Output("analysis_report","children", allow_duplicate=True),
    Input("analysis-show-editor","value"),
    prevent_initial_call=True
)
def activate_report(val):
        
    if val != None:

        l = markdown_to_dash(config.report)

        if val:

            d = [
                    dbc.Row([
                            dbc.Col(id="analysis-markdown", children=l),
                            dbc.Col(
                                dash_editor_components.PythonEditor(
                                    id='analysis-editor',
                                    value=config.report
                                )
                            ),
                    ])
                ]
            
        else:

            d = [
                    dbc.Col(id="analysis-markdown", children=l),
                ]

        return d

    else:

        raise PreventUpdate() 

@app.callback(
    Output("analysis-markdown","children", allow_duplicate=True),
    Input("analysis-editor","value"),
    prevent_initial_call=True
)
def editor(value):

    if value:
        config.report = value

        l = markdown_to_dash(value)

        return l
    
    else:

        raise PreventUpdate()
    
@app.callback(
    Output("dumb","children", allow_duplicate=True),
    Input("analysis-save-report","n_click"),
    prevent_initial_call=True
)
def editor(value):

    if value:

        save_graph()

        return ""
    
    else:

        raise PreventUpdate()
    
@app.callback(
    Output("analysis_info_plot","children", allow_duplicate=True),
    Input("analysis_plot_dropdown","value"),
    prevent_initial_call=True
)
def plot_info(value):

    if value != None:

        print("holi")

        docs = config.methods_plot[value]["docs"]

        print(docs)

        return html.Pre(docs, style={"white-space":"pre-wrap"})
    
    else:

        raise PreventUpdate()
