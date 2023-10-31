from dash import dcc, html, dash_table, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto

import os

from app import app

from general import *

for i in os.listdir("./methods"):
    if i.endswith(".py"):
        f = f"from methods.{i[:-3]} import *"
        exec(f)

############################################################################################################################################
############################################################################################################################################
# Layout
############################################################################################################################################
############################################################################################################################################
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
                    dbc.ModalHeader("New",close_button=False),
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
                        dbc.Button(
                            children="Delete",
                            id="analysis_delete_button", 
                            style={"background-color":"red"}, 
                        ),
                    ),
                    dbc.Col(
                        dcc.Dropdown(
                            id='graph_dropdown_analysis',
                            value=[i for i in config.methods.keys() if i != "Raw"][0],
                            options=[i for i in config.methods.keys() if i != "Raw"],
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
                cyto.Cytoscape(
                    id='graph_cytoscape',
                    layout={'name': 'preset'},
                    style={'width': '100%', 'height': '500px'},
                    elements=config.graph,
                    userZoomingEnabled=True,  # Disable zooming
                    userPanningEnabled=False,  # Disable panning
                    stylesheet=GRAPHSTYLESHEET,
                ),
                style={'width': '100wh', 'height':'120wh', 'margin': 'auto'}, #vh/wh
            ),
            dbc.Row(id = 'analysis_args',
                    children = [],
                    # width='50%',
                    style=ARGUMENTBOXSTYLE
                    ),
            dbc.Row(id='analysis_postargs',
                    children=[],  
                    # width='50%',
                    style=ARGUMENTBOXSTYLE
                    ),
            dbc.Row(id='analysis_plotargs',
                    children=[],  
                    # width='50%',
                    style=ARGUMENTBOXSTYLE
                    ),
            dbc.Row(id='analysis_plot',
                    children=[]
                    )
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
                elif populate == "postexecution":
                    pos = get_node_pos(config.selected)
                    config.graph[pos]["data"]["parameters"][args[i]["name"]] = value
                elif populate == "plot":
                    pos = get_node_pos(config.selected)
                    config.graph[pos]["data"]["plot"][args[i]["name"]] = value

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

def method_create_pars(args_list):

    method_args = config.methods[get_node(config.selected)["data"]["method"]]["args"]
    for i in args_list:
        args = parameters_eval(method_args[i], i)
        if i == "execution":
            pos = get_node_pos(config.selected)

def make_arguments(id, arg_list, loaded_args={}, add_execution_button=True, add_header="args"):

    #Avoid updates after creating the inputs
    arglist = parameters_eval(arg_list)
    for i in arglist:
        config.block_callback[i["name"]] = True

    l = [
        dbc.Col(
            dbc.Row(html.H1("Additional arguments")),
        ),
    ]
    if add_header == "args":
        l = [
            dbc.Row(
                [
                    dbc.Col(
                        html.H1(config.selected, id="analysis_name")
                    ),
                ],
                align="center"
            )
        ]
    elif add_header == "plot":
        l = [
            dbc.Col(
                dbc.Row(html.H1("Plot arguments")),
            ),
            dbc.Col(
                dbc.Row(
                    dbc.Button("Save Image",id="analysis_saveimage_button", class_name="btn btn-primary btn-sm"),
                )
            )
        ]

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
                    None

                input = dbc.Input(
                            id="analysis_"+str(arg["name"]),
                            **arg["properties"]
                        )

            elif arg["input"] == "Dropdown":

                try:
                    arg["properties"]["value"] = loaded_args[arg["name"]]
                except:
                    None

                input = dcc.Dropdown(
                            id="analysis_"+str(arg["name"]),
                            **arg["properties"]
                        )
                
            elif arg["input"] == "BooleanSwitch":

                try:
                    arg["properties"]["on"] = loaded_args[arg["name"]]
                except:
                    None

                input = daq.BooleanSwitch(
                            id="analysis_"+str(arg["name"]),
                            **arg["properties"]
                        )

            elif arg["input"] == "AgTable":

                try:
                    arg["properties"]["data"] = loaded_args[arg["name"]]
                except:
                    None

                if "deleteRows" in arg.keys():
                    d = arg["deleteRows"]
                else:
                    d = False

                input = [
                    dbc.Row(
                        html.Div(
                                children= ag_table("analysis_"+str(arg["name"]), arg["properties"]["header"], arg["properties"]["data"], deleteRows=d),
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

def load_node(name):

    if name == "Raw":
        return [], [], [], []

    method = get_node(name)['data']['method']
    args = get_node(name)['data']['parameters']
    plot_args = get_node(name)['data']['plot']
    l = make_arguments(method, config.methods[method]["args"]["execution"], args)
    l2 = []
    l3 = []
    p = []
    if get_node(name)['data']['computed']:
        l2 = make_arguments(method, config.methods[method]["args"]["postexecution"], args, add_execution_button=False, add_header="postargs")
        l3 = make_arguments(method, config.methods[method]["args"]["plot"], plot_args, add_execution_button=False, add_header="plot")
        p = config.methods[method]["plot"]()
            
    return l, l2, l3, p

############################################################################################################################################
############################################################################################################################################
# Callbacks
############################################################################################################################################
############################################################################################################################################

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
                input = f"Input('analysis_{i['name']}','value'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = f"config.graph[get_node_pos(config.selected)]['data']['parameters']['{i['name']}'] = data"
                up_plot = f"config.graph[get_node_pos(config.selected)]['data']['plot']['{i['name']}'] = data"
                args = "data"
            elif i['input'] == 'Dropdown':
                input = f"Input('analysis_{i['name']}','value'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = f"config.graph[get_node_pos(config.selected)]['data']['parameters']['{i['name']}'] = data"
                up_plot = f"config.graph[get_node_pos(config.selected)]['data']['plot']['{i['name']}'] = data"
                args = "data"
            elif i['input'] == 'BooleanSwitch':
                input = f"Input('analysis_{i['name']}','on'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = f"config.graph[get_node_pos(config.selected)]['data']['parameters']['{i['name']}'] = data"
                up_plot = f"config.graph[get_node_pos(config.selected)]['data']['plot']['{i['name']}'] = data"
                args = "data"
            elif i['input'] == 'AgTable':
                input = f"Input('analysis_{i['name']}','cellValueChanged'), State('analysis_{i['name']}','rowData'),"
                up = f"if cell != None: data[cell['rowIndex']] = cell['data']; config.active_node_parameters['{i['name']}'] = data"
                up_post = f"if cell != None: data[cell['rowIndex']] = cell['data']; config.graph[get_node_pos(config.selected)]['data']['parameters']['{i['name']}'] = data"
                up_plot = f"if cell != None: data[cell['rowIndex']] = cell['data']; config.graph[get_node_pos(config.selected)]['data']['plot']['{i['name']}'] = data"
                args = "cell, data"

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
    if '{i['name']}' not in [i['name'] for i in config.methods[method]["args"]["execution"]]:
        raise PreventUpdate()

    if config.block_callback['{i['name']}']:
        config.block_callback['{i['name']}'] = False
        raise PreventUpdate()

    {up}
    clean_arguments('{i['name']}', config.methods[method]["args"]["execution"], config.active_node_parameters)
    args_object = make_arguments(method, config.methods[method]["args"]["execution"], config.active_node_parameters)
    config.block_callback['{i['name']}'] = False

    return args_object
"""

            exec(add_function, globals(), locals())

#Postexecution and Plot
            add_function = f"""
@app.callback(
    Output("analysis_postargs","children", allow_duplicate=True),
    Output("analysis_plotargs","children", allow_duplicate=True),
    Output("analysis_plot","children", allow_duplicate=True),
    {input}
    prevent_initial_call=True
)
def change_postparameter_{i['name']}({args}):

    method = get_node(config.selected)["data"]["method"]
    if '{i['name']}' in [i['name'] for i in config.methods[method]["args"]["postexecution"]]:

        if '{i['name']}' in config.block_callback.keys():
            if config.block_callback['{i['name']}']:
                config.block_callback['{i['name']}'] = False
                raise PreventUpdate()

        {up_post}
        set_parameters(dict({i['name']}=data), 'parameters')
    
    elif '{i['name']}' in [i['name'] for i in config.methods[method]["args"]["plot"]]:
    
        if '{i['name']}' in config.block_callback.keys():
            if config.block_callback['{i['name']}']:
                config.block_callback['{i['name']}'] = False
                raise PreventUpdate()

        {up_plot}
        set_parameters(dict({i['name']}=data), 'parameters')

    else:
    
        raise PreventUpdate()

    method = get_node(config.selected)["data"]["method"]

    post_args = deepcopy(config.methods[method]["args"]["postexecution"])
    deactivate = False
    for i in post_args:
        config.block_callback[i['name']] = True
    config.block_callback['{i['name']}'] = False
    clean_arguments('{i['name']}', config.methods[method]["args"]["postexecution"], get_node(config.selected)['data']['parameters'])
    post_args_object = make_arguments(method, post_args, loaded_args=get_node(config.selected)['data']['parameters'], add_execution_button=False, add_header="postargs")

    plot_args = deepcopy(config.methods[method]["args"]["plot"])
    deactivate = False
    clean_arguments('{i['name']}', config.methods[method]["args"]["plot"], get_node(config.selected)['data']['plot'])
    plot_args_object = make_arguments(method, plot_args, loaded_args=get_node(config.selected)['data']['plot'], add_execution_button=False, add_header="plot")
    config.block_callback['{i['name']}'] = False

    plot = config.methods[method]['plot']()

    return post_args_object, plot_args_object, plot
"""

            exec(add_function, globals(), locals())

#Delete row ag
            if i['input'] == "AgTable" and 'deleteRows' in i.keys():

                if i['deleteRows']:

                    add_function = f"""
@app.callback(
    Output('analysis_{i['name']}', "rowData", allow_duplicate=True),
    Input('analysis_{i['name']}', "virtualRowData"),
    State('analysis_{i['name']}', "rowData"),
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
    Output("analysis_{i['name']}","rowData", allow_duplicate=True),
    Input("analysis_{i['name']}_button","n_clicks"),
    State("analysis_{i['name']}","rowData"),
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

@app.callback(
    Output('new-modal', 'is_open', allow_duplicate=True),
    Output('new-dropdown-before', 'value', allow_duplicate=True),
    Output('new-dropdown-before', 'options', allow_duplicate=True),
    Input('graph_new_button', 'n_clicks'),
    prevent_initial_call=True
)
def graph_new_modal(_):

    return True, config.selected, node_names()

@app.callback(
    Output('new-modal', 'is_open', allow_duplicate=True),
    Input('new-cancel', 'n_clicks'),
    prevent_initial_call=True
)
def graph_new_node(_):

    return False

@app.callback(
    Output('new-modal', 'is_open', allow_duplicate=True),
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('graph_table', 'data', allow_duplicate=True),
    Output('analysis_args', 'children', allow_duplicate=True),
    Output('analysis_postargs', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    [
        Input('new-proceed', 'n_clicks')
    ],
    [
        State('new-dropdown-before', 'value'),
        State('new-dropdown-after', 'value'),
        State('graph_dropdown_analysis', 'value'),
        State('graph_cytoscape', 'layout'),
    ],
    prevent_initial_call=True
)
def graph_new_node(_, input, output, method, cytoscape):

    #Make new name
    name = method
    nodes = node_names()
    count = 0
    while name in nodes:
        count += 1
        name = f"{method}_{str(count)}"

    node = {
        'data': {
            'id': name, 
            'name': method,
            'method':method, 
            'type':config.methods[method]["properties"]['type'], 
            'color': "blue",
            'h5ad_file':None,
            'image': '',
            'computed':False,
            'opacity':.3,
            'summary':'', 
            'parameters':{"input":input},
            'plot':{}
        }, 
        'position':{'x':config.max_x + 30,'y':0},
    }
    #make active node that is the one that will be presented
    # nodes_update_pos(config.graph,cytoscape)

    if input != config.h5ad_file: #Load appropiate node
        # save_adata()
        innode = get_node(input)["data"]
        load_adata(f"{config.analysis_folder}/h5ad/{innode['h5ad_file']}")

    #Create node
    config.graph.append(node)
    #Create uns
    config.adata.uns[name] = {}
    config.adata.uns[name]["plot"] = {}
    #Create block callbacks
    config.block_callback = {}
    #Selected
    unselect_node(config.selected)
    config.selected = name
    select_node(config.selected)
    #Create parameters
    method_create_pars(["execution"])
    #Create active
    config.active_node_parameters["input"] = input
    set_parameters(config.active_node_parameters, "parameters")
    # set_active_node_parameters({i["name"]:get_value(i) for i in args.copy()})
    #Assign new file
    innode = get_node(input)["data"]
    pos = get_node_pos(config.selected)
    if config.methods[innode["method"]]["properties"]["make_new_h5ad"]:
        config.graph[pos]["data"]["h5ad_file"] = new_h5ad_file()
    else:
        config.graph[pos]["data"]["h5ad_file"] = innode["h5ad_file"]

    edge_add(input, config.selected)

    if output:

        deactivate_downstream(output)
        deactivate_node(output)
        node_reassign_input(output,name)

    #load interactive plots
    l, l2, l3, p = load_node(name)

    return False, config.graph, graph2table(), l, l2, l3, p 

#Execute analysis button
@app.callback(
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('execute-computed-modal', 'is_open', allow_duplicate=True),
    Output('analysis_postargs', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    [
        Input('analysis_execute_button', 'n_clicks')
    ],
    State('execute-computed-modal', 'is_open'),
    State('analysis_plot', 'children'),
    prevent_initial_call=True
)
def execute(n_clicks, warning_computed, plot):

    if n_clicks != None:

        for node in get_node_ancestors(config.selected):
            
            if not node['data']['computed']:
            
                return config.graph, True, [], [], []

        #Get incoming input, old input and method
        input = config.active_node_parameters["input"]

        method = get_node(config.selected)["data"]["method"]

        #Get erguments of incoming
        inputArgs = get_args(input)

        #Execute code
        innode = get_node(input)["data"]
        pos = get_node_pos(config.selected)
        if innode["h5ad_file"] != config.h5ad_file: #Load appropiate node
            load_adata(f"{config.analysis_folder}/h5ad/{innode['h5ad_file']}")

        # if config.methods[innode["method"]]["properties"]["make_new_h5ad"]:
        #     config.graph[pos]["data"]["h5ad_file"] = new_h5ad_file()
        #     save = True
        # else:
        #     config.graph[pos]["data"]["h5ad_file"] = innode["h5ad_file"]
        #     save = False

        config.h5ad_file = config.graph[pos]["data"]["h5ad_file"]      
        config.methods[method]["function"](config.adata, config.active_node_parameters)

        #save parameters
        set_parameters(config.active_node_parameters, 'parameters')
        method_create_pars(["postexecution","plot"])

        #Activate node
        activate_node(config.selected)
        deactivate_downstream(config.selected)

        node_data = get_node(config.selected)['data']
        post_args_object = make_arguments(node_data['method'], config.methods[node_data['method']]["args"]["postexecution"], add_execution_button=False, add_header="postargs")
        plot_args_object = make_arguments(node_data['method'], config.methods[node_data['method']]["args"]["plot"], add_execution_button=False, add_header="plot")
        plot = config.methods[node_data['method']]["plot"]()

        save_adata()
        save_graph()

    else:

        raise PreventUpdate()

    return config.graph, warning_computed, post_args_object, plot_args_object, plot

#Delete analysis button
@app.callback(
    Output('delete-modal', 'is_open'),
    Output('delete-message', 'children'),
    Input('analysis_delete_button', 'n_clicks'),
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
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('analysis_args', 'children', allow_duplicate=True),
    Output('analysis_postargs', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    Output('delete-modal', 'is_open', allow_duplicate=True),
    Input('delete-proceed', 'n_clicks'),
    State('graph_dropdown_load', 'value'),
    prevent_initial_call=True
)
def delete_confirmation(n_clicks, val):
    
    modal = True
    if n_clicks != None or  val != None:

        if val != "Raw":

            deactivate_downstream(val)
            node_rm(val)

            if val == config.selected:
                config.selected = 'Raw'
            
            modal = False

            l, l2, l3, p = load_node(config.selected)

            return config.graph, l, l2, l3, p, modal
        
        else:

            l, l2, l3, p = load_node(config.selected)
            modal = False

            return config.graph, l, l2, l3, p, modal
    
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
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('analysis_args', 'children', allow_duplicate=True),
    Output('analysis_postargs', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
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
        if node["h5ad_file"] != config.h5ad_file:
            #reload file previous to the model
            load_adata(f"{config.analysis_folder}/h5ad/{node['h5ad_file']}")
            config.h5ad_file = node["h5ad_file"]
            # adapt_adata_loaded()

        #Create block callbacks
        config.block_callback = {}
        #Create active
        pos = get_node_pos(name)
        config.active_node_parameters = config.graph[pos]["data"]["parameters"]
        #Selected
        unselect_node(config.selected)
        config.selected = name
        select_node(config.selected)

        # #Create parameters
        # method_create_pars(["execution"])
        # # set_active_node_parameters({i["name"]:get_value(i) for i in args.copy()})

        l, l2, l3, p = load_node(name)

        return config.graph, l, l2, l3, p
    
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
    Output("dumb","children",allow_duplicate=True),
    Input("analysis_saveimage_button","n_clicks"),
    State("analysis_plot","children"),
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
