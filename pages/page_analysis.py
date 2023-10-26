from dash import dcc, html, dash_table, Input, Output, State
from dash.exceptions import PreventUpdate
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
import dash_ag_grid as dag

import os
import json

from app import app

from . import config
from .functions import *

for i in os.listdir("./pages/methods"):
    if i.endswith(".py"):
        f = f"from .methods.{i[:-3]} import *"
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
def save_graph():

    unselect_node(config.selected)
    file = f"{config.analysis_folder}/analysis.json"
    with open(file,"w") as outfile:
        json_object = json.dumps(config.graph, indent=4, cls=NpEncoder)
        outfile.write(json_object)
    select_node(config.selected)

def save_adata():

    pos = get_node_pos(config.selected)
    file = f"{config.analysis_folder}/h5ad/{config.graph[pos]['data']['h5ad_file']}"
    config.adata.write(file)

def save_report():
    
    file = f"{config.analysis_folder}/report/report.md"
    with open(file,"w") as outfile:
        outfile.write(config.report)

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
    config.graph[pos]['data']['parameters'] = deepcopy(config.active_node_parameters)

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

def unselect_node(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'.png'
    
def select_node(name):
    pos = get_node_pos(name)
    node = get_node(name)

    config.graph[pos]['data']['image'] = '../assets/'+node['data']['type']+'_selected.png'

def list_observables():
    l = list(config.adata.obs.columns.values)
    ancestors = get_nodes()
    for i in ancestors:
        if "obs" in i["data"].keys():
            for j in i["data"]["obs"].keys():
                l += [f"{i['data']['name']}--{j}"]

    return l

def prevent_race(name,computed=True,method=True):

    node = get_node(config.selected)
    if not node['data']['computed'] and computed:
        raise PreventUpdate()

    if node['data']['method'] != name and method:
        raise PreventUpdate()
    
def is_computed():

    return get_node(config.selected)["data"]["computed"]

def graph2table():
    return [{"Name":i['data']['id'],"Type":i['data']['type'],"Method":i['data']['method']} for i in get_nodes() if i['data']['id'] != 'Raw']

# def set_active_node_parameters(args):
#     config.active_node_parameters = args

def set_parameters(args, arg_type):

    pos = get_node_pos(config.selected)
    for i,j in args.items():
        config.graph[pos]["data"][arg_type][i] = j

    for i,j in args.items():
        config.adata.uns[config.selected][arg_type][i] = j
    config.adata.uns[config.selected]["scinteractive"] = True
    config.adata.uns[config.selected]["method"] = get_node(config.selected)["data"]["method"]

for i in os.listdir("./pages/methods"):
    if i.endswith(".py"):
        f = f"from .methods.{i[:-3]} import *"
        exec(f)

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

            if count == 0 and\
               populate == "execution":
                config.active_node_parameters[args[i]["name"]] = args[i]["value"]
                # config.adata.uns[config.selected]["parameters"][args[i]["name"]] = args[i]["value"]
                # pos = get_node_pos(config.selected)
                # config.graph[pos]["data"]["parameters"][args[i]["name"]] = args[i]["value"]
            elif count == 0 and\
               populate == "postexecution":
                config.adata.uns[config.selected]["parameters"][args[i]["name"]] = args[i]["value"]
                pos = get_node_pos(config.selected)
                config.graph[pos]["data"]["parameters"][args[i]["name"]] = args[i]["value"]
            elif count == 0 and\
               populate == "plot":
                config.adata.uns[config.selected]["plot"][args[i]["name"]] = args[i]["value"]
                pos = get_node_pos(config.selected)
                config.graph[pos]["data"]["plot"][args[i]["name"]] = args[i]["value"]

    elif type(args) == dict:
        if "function" in args.keys():
            return fvalue(args)
        else:
            for i,j in args.items():
                args[i] = parameters_eval(j, populate, count+1)

    return args

def parameters2args(args, populate):

    args = parameters_eval(args, populate)

    return {i["name"]:i["value"] for i in args}    

def method_create_pars(args_list):

    method_args = config.methods[get_node(config.selected)["data"]["method"]]["args"]
    for i in args_list:
        parameters_eval(method_args[i], i)

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
        config.adata.obsm[config.selected] = args["obsm"]
            
    #uns
    if "uns" in args.keys():
        config.adata.uns[config.selected] = args["uns"]

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

        arg = arg.copy()

        if fvisible(arg):

            for i,j in arg.items():
                arg[i] = fvalue(j)

            if loaded_args != {}:
                try:
                    value = loaded_args[arg["name"]]
                except:
                    value = None
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

def new_h5ad_file():

    count = 1
    file_name = f"Analysis_{count}.h5ad"
    while os.path.exists(f"{config.analysis_folder}/h5ad/{file_name}"):
        count += 1
        file_name = f"Analysis_{count}.h5ad"

    return file_name

def clean_h5ad():

    l = os.listdir(f"{config.analysis_folder}/h5ad")
    f = [i["data"]["h5ad_file"] for i in get_nodes()]
    for file in l:
        if file not in f:
            os.remove(f"{config.analysis_folder}/h5ad/{file}")
        
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
                up_post = ""
                up_plot = ""
                args = "data"
            elif i['input'] == 'Dropdown':
                input = f"Input('analysis_{i['name']}','value'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = ""
                up_plot = ""
                args = "data"
            elif i['input'] == 'BooleanSwitch':
                input = f"Input('analysis_{i['name']}','on'),"
                up = f"config.active_node_parameters['{i['name']}'] = data"
                up_post = ""
                up_plot = ""
                args = "data"
            elif i['input'] == 'AgTable':
                input = f"Input('analysis_{i['name']}','cellValueChanged'), State('analysis_{i['name']}','rowData'),"
                up = f"if cell != None: data[cell['rowIndex']] = cell['data']; config.active_node_parameters['{i['name']}'] = data"
                up_post = "if cell != None: data[cell['rowIndex']] = cell['data']"
                up_plot = up_post
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
    clean_arguments('{i['name']}', config.methods[method]["args"]["postexecution"], config.adata.uns[config.selected]['parameters'])
    post_args_object = make_arguments(method, post_args, loaded_args=config.adata.uns[config.selected]['parameters'], add_execution_button=False, add_header="postargs")

    plot_args = deepcopy(config.methods[method]["args"]["plot"])
    deactivate = False
    for i in plot_args:
        config.block_callback[i['name']] = True
    config.block_callback['{i['name']}'] = False
    clean_arguments('{i['name']}', config.methods[method]["args"]["plot"], config.adata.uns[config.selected]['plot'])
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
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('graph_table', 'data', allow_duplicate=True),
    Output('analysis_args', 'children', allow_duplicate=True),
    Output('analysis_postargs', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    [
        Input('graph_new_button', 'n_clicks')
    ],
    [
        State('graph_dropdown_analysis', 'value'),
        State('graph_cytoscape', 'layout'),
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
            'parameters':{},
            'plot':{}
        }, 
        'position':{'x':config.max_x + 30,'y':0},
    }
    #make active node that is the one that will be presented
    # nodes_update_pos(config.graph,cytoscape)

    #Create node
    config.graph.append(node)
    #Create uns
    config.adata.uns[name] = {}
    config.adata.uns[name]["parameters"] = {}
    config.adata.uns[name]["plot"] = {}
    #Create block callbacks
    config.block_callback = {}
    #Create active
    config.active_node_parameters = {}
    #Selected
    unselect_node(config.selected)
    config.selected = name
    select_node(config.selected)
    #Create parameters
    method_create_pars(["execution"])
    # set_active_node_parameters({i["name"]:i["value"] for i in args.copy()})

    #load interactive plots
    l, l2, l3, p = load_node(name)

    return config.graph, graph2table(), l, l2, l3, p 

#Execute analysis button
@app.callback(
    Output('graph_cytoscape', 'elements', allow_duplicate=True),
    Output('execute-input-modal', 'is_open', allow_duplicate=True),
    Output('execute-computed-modal', 'is_open', allow_duplicate=True),
    Output('analysis_postargs', 'children', allow_duplicate=True),
    Output('analysis_plotargs', 'children', allow_duplicate=True),
    Output('analysis_plot', 'children', allow_duplicate=True),
    [
        Input('analysis_execute_button', 'n_clicks')
    ],
    State('execute-input-modal', 'is_open'),
    State('execute-computed-modal', 'is_open'),
    State('analysis_plot', 'children'),
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
                try: 
                    old_input = get_node(config.selected)["data"]["parameters"]["input"]
                except:
                    old_input = "Raw"
                method = get_node(config.selected)["data"]["method"]

                #Get erguments of incoming
                inputArgs = get_args(input)

                #Execute code
                innode = get_node(input)["data"]
                pos = get_node_pos(config.selected)

                if innode["h5ad_file"] != config.h5ad_file: #Load appropiate node
                    del config.adata.uns[config.selected]
                    config.adata = sc.read(f"{config.analysis_folder}/h5ad/{innode['h5ad_file']}")
                    config.adata.uns[config.selected] = {"parameters":deepcopy(config.active_node_parameters),"plot":{}}

                if config.methods[innode["method"]]["properties"]["make_new_h5ad"]:
                        
                    clean_h5ad()
                    config.graph[pos]["data"]["h5ad_file"] = new_h5ad_file()
                    save = True

                else:

                    config.graph[pos]["data"]["h5ad_file"] = innode["h5ad_file"]
                    save = False

                config.h5ad_file = config.graph[pos]["data"]["h5ad_file"]      

                outputArgs = config.methods[method]["function"](config.adata, inputArgs, config.active_node_parameters)

                # Remove old edge if any and add new
                edge_rm(old_input, config.selected)
                edge_add(input, config.selected)

                #save parameters
                set_output(outputArgs)
                set_parameters(config.active_node_parameters, 'parameters')

                method_create_pars(["postexecution","plot"])

                if save: #Save if needed to generate a new file
                    adapt_adata_saving()
                    save_adata()
                    adapt_adata_loaded()

                #Activate node
                activate_node(config.selected)
                deactivate_downstream(config.selected)

        post_args_object = {}
        plot = []
        node_data = get_node(config.selected)['data']
        if node_data['computed']:

            post_args_object = make_arguments(node_data['method'], config.methods[node_data['method']]["args"]["postexecution"], add_execution_button=False, add_header="postargs")
            plot_args_object = make_arguments(node_data['method'], config.methods[node_data['method']]["args"]["plot"], add_execution_button=False, add_header="plot")

            plot = config.methods[node_data['method']]["plot"]()

    else:

        raise PreventUpdate()

    return config.graph, warning_input, warning_computed, post_args_object, plot_args_object, plot

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
    prevent_initial_call=True
)
def delete_confirmation(n_clicks):
    
    modal = True
    if n_clicks != None:

        name = config.selected
        deactivate_downstream(name)
        node_rm(name)
        clean_h5ad()

        config.selected = 'Raw'
        modal = False

        l, l2, l3, p = load_node(config.selected)

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

        node = get_node(name)["data"]
        if node["h5ad_file"] != config.h5ad_file:
            #reload file previous to the model
            config.adata = sc.read(f"{config.analysis_folder}/h5ad/{node['h5ad_file']}")
            config.h5ad_file = node["h5ad_file"]
            adapt_adata_loaded()

        #Create block callbacks
        config.block_callback = {}
        #Create active
        config.active_node_parameters = {}
        #Selected
        unselect_node(config.selected)
        config.selected = name
        select_node(config.selected)

        # #Create parameters
        # method_create_pars(["execution"])
        # # set_active_node_parameters({i["name"]:i["value"] for i in args.copy()})

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

        save_graph()

        adapt_adata_saving()
        save_adata()
        adapt_adata_loaded()

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
