import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import dash
import scrublet
from scipy.stats import mode
import dash_ag_grid as dag
import base64
import datetime
import io

from general import *

def args_marker_genes():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_layer = [i for i in config.adata.layers.keys()]

    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Use the indicated neighbors.",
            "value":None,
            "clearable":False,
            "options":options
        },
    ]

def f_marker_genes(name_analysis, kwargs):

    pos = get_node_pos(name_analysis)
    config.graph[pos]["data"]["plotting"] = {"markers":{"Types":{"0":"Example TypeX"},"Genes":{"0":"[ExampleGene1, ExampleGene2]"}}, "var":config.adata.var.columns[0]}

def rm_marker_genes(name_analysis):

    return

def rename_marker_genes(name_analysis, name_new_analysis):

    return

def plot_marker_genes(name_analysis):
    
    node = get_node(name_analysis)

    if not node["data"]["computed"]:
        return []

    m = node["data"]["plotting"]["markers"]["Genes"]
    genes = []
    for i,j in m.items():
        genes += j[1:-1].split(",")
    
    name = node["data"]["parameters"]["input"]
    # clusters = config.adata.uns[name]["scores"].dtype.names

    data_array, xorder, yorder = de_markers2array(config.adata.uns[name],genes,node["data"]["plotting"]["var"])

    table = plot_table(pd.DataFrame(node["data"]["plotting"]["markers"]),False,False,False,editable=True, rowdeletable=True)
    plot = plot_scattermap(data_array)

    plot = [
            table,
            html.Div(
            dcc.Upload(
                id = "upload_file",
                children=html.Div([
                    'Drag and Drop or Select Gene Files'
                ]),
                style={
                    'width': '100%',
                    'height': '60px',
                    'lineHeight': '60px',
                    'borderWidth': '1px',
                    'borderStyle': 'dashed',
                    'borderRadius': '5px',
                    'textAlign': 'center',
                    'margin': '10px'
                },
                multiple=False
            ),
            ),
            dcc.Dropdown(
                id = "var_name",
                value = node["data"]["plotting"]["var"],
                options = config.adata.var.columns.values
            ),
            dcc.Graph(figure={
                            "data":[plot],
                            "layout":{
                                "xaxis":{
                                    "tickmode": "array",
                                    "tickvals": np.arange(0,len(xorder),1),
                                    "ticktext": xorder
                                },
                                "yaxis":{
                                    "tickmode": "array",
                                    "tickvals": np.arange(0,len(yorder),1),
                                    "ticktext": yorder
                                }
                            }
                    })
    ]

    return plot

# @app.callback(dash.Output('analysis_plot', 'children'),
#           dash.Input('upload_file', 'contents'),
#           dash.Input('var_name', 'value'),
#           dash.State('upload_file', 'filename'),
#           dash.State('upload_file', 'last_modified'))
# def update_output(list_of_contents, var, list_of_names, list_of_dates):

#     m = html.Div()

#     plot = [
#             m,
#             html.Div(
#             dcc.Upload(
#                 id = "upload_file",
#                 children=html.Div([
#                     'Drag and Drop or Select Files'
#                 ]),
#                 style={
#                     'width': '100%',
#                     'height': '60px',
#                     'lineHeight': '60px',
#                     'borderWidth': '1px',
#                     'borderStyle': 'dashed',
#                     'borderRadius': '5px',
#                     'textAlign': 'center',
#                     'margin': '10px'
#                 },
#                 multiple=False
#             ),
#             )
#     ]

#     pos = get_node_pos(config.selected)
#     config.graph[pos]["data"]["plotting"]["var"] = var

#     if list_of_contents != None:

#         content_type, content_string = list_of_contents.split(',')

#         decoded = base64.b64decode(content_string)
#         try:
#             if 'csv' in list_of_names:
#                 # Assume that the user uploaded a CSV file
#                 df = pd.read_csv(
#                     io.StringIO(decoded.decode('utf-8')), sep=":",header=None)
                
#                 # df = df.astype(str)
#                 df.columns = ["Type","Genes"]
    
#                 config.graph[pos]["data"]["plotting"]["markers"] = df.to_dict()

#                 plot = plot_marker_genes(config.selected)

#         except Exception as e:
#             None

#     elif len(get_node(config.selected)["data"]["plotting"]["markers"]) != 0:

#         plot = plot_marker_genes(config.selected)
    
#     return plot

# @app.callback(
#     dash.Output('analysis_plot', 'children', allow_duplicate=True),
#     dash.Input("table", "cellRendererData"),
#     prevent_initial_call = True
# )
# def delete(n):

#     prevent_race("marker_genes")

#     node = get_node(config.selected)
#     pos = get_node_pos(config.selected)

#     if n:
#         if n["rowIndex"] < len(node["data"]["plotting"]["markers"]):
#             del config.graph[pos]["data"]["plotting"]["markers"]["Type"][n["rowId"]]
#             del config.graph[pos]["data"]["plotting"]["markers"]["Genes"][n["rowId"]]

#     l = np.arange(0,len(config.graph[pos]["data"]["plotting"]["markers"]["Type"]),1)
#     k = [int(i) for i in config.graph[pos]["data"]["plotting"]["markers"]["Type"].keys()]
#     l_1 = [i for i in l if i not in k]
#     k_1 = [i for i in k if i not in l]
#     for i,j in zip(l_1, k_1):
#         config.graph[pos]["data"]["plotting"]["markers"]["Type"][str(i)] = config.graph[pos]["data"]["plotting"]["markers"]["Type"][str(j)]
#         del config.graph[pos]["data"]["plotting"]["markers"]["Type"][str(j)]

#         config.graph[pos]["data"]["plotting"]["markers"]["Genes"][str(i)] = config.graph[pos]["data"]["plotting"]["markers"]["Genes"][str(j)]
#         del config.graph[pos]["data"]["plotting"]["markers"]["Genes"][str(j)]

#     plot = plot_marker_genes(config.selected)

#     return plot

# @app.callback(
#     dash.Output('analysis_plot', 'children', allow_duplicate=True),
#     dash.Input("table", "cellValueChanged"),
#     prevent_initial_call = True
# )
# def delete(n):

#     prevent_race("marker_genes")

#     node = get_node(config.selected)
#     pos = get_node_pos(config.selected)

#     if n != None:

#         if n["rowIndex"] < len(node["data"]["plotting"]["markers"]):
#             config.graph[pos]["data"]["plotting"]["markers"]["Type"][n["rowId"]] = n["data"]["Type"]
#             config.graph[pos]["data"]["plotting"]["markers"]["Genes"][n["rowId"]] = n["data"]["Genes"]
#         else:
#             config.graph[pos]["data"]["plotting"]["markers"]["Type"][str(n["rowIndex"])] = n["data"]["Type"]
#             config.graph[pos]["data"]["plotting"]["markers"]["Genes"][str(n["rowIndex"])] = n["data"]["Genes"]

#     else:
    
#         raise PreventUpdate()

#     plot = plot_marker_genes(config.selected)

#     return plot