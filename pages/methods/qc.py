import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash
from scipy.stats import mode

from ..arguments import *
from ..functions import *
from ..graph import *
from ..plots import *

from app import app

from .. import config

args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"QCTable",
            "name":"measure",
            "description":"Metrics of quality control to compute",
            "value":[
                {
                    "name":"counts",
                    "var":" ",
                    "pattern":" ",
                    "style":"counts",
                    "genes":"all",
                },
                {
                    "name":"n_expressed_genes",
                    "var":" ",
                    "pattern":" ",
                    "style":"n_expressed_genes",
                    "genes":"all",
                }
            ],
        },
    ],

    "postexecution" : [
        ARGBATCH,
        {
            "input":"Dropdown",
            "name":"plot_style",
            "description":"Chose between violin or 2D scatterplot.",
            "value":"violin",
            "clearable":False,
            "options":["violin","scatter"],
        },
        {
            "input":"Dropdown",
            "name":"plot_measure",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":{"function": "qc_measures()[0]"},
            "clearable":False,
            "options":{"function": "qc_measures()"},
            "visible":{"function": "config.adata.uns[config.selected]['parameters']['plot_style'] == 'scatter'"}
        },
        {
            "input":"AgTable",
            "name":"thresholds",
            "description":"Thresholds to apply to the data",
            "header":[
                { "headerName": "Metric", "field":"metric", "editable": False },
                { "headerName": "Cluster", "field":"cluster", "editable": False },
                { "headerName": "Min", "field":"min", "editable": True },
                { "headerName": "Max", "field":"max", "editable": True },
                { "headerName": "Plot resolution", "field":"resolution", "editable": True },
            ],
            "value":{"function":"qc_data()"},
        },
    ]

}

def qc_measures():
    selected = config.selected
    adata = config.adata

    return [i["name"] for i in adata.uns[selected]["parameters"]["measure"]]

def qc_data():

    selected = config.selected
    adata = config.adata

    a = adata.uns[selected]["parameters"]

    # #Reset table
    # if "batch" in modified_arg.keys():
    #     if "batch" in config.adata.uns[config.selected]["parameters"].keys():
    #         del config.adata.uns[config.selected]["parameters"]["batch"]

    measures = [i["name"] for i in adata.uns[selected]["parameters"]["measure"]]

    measure = a["measure"][0]
    if "plot_measure" in a.keys():
        measure = a["plot_measure"]

    data = []
    for i in a["measure"]:
        if i["name"] in measures:
            i = i["name"]
            j = f"{config.selected}--{measure}"
            data.append({"metric":i, "cluster":None,"min":adata.obs[j].min(),"max":adata.obs[j].max(),"resolution":adata.obs[j].std()})

    return data

def f_qc(adata, inputArgs, kwargs):
        
    X = inputArgs["obsm"]

    d = {"obs":{}}
    for l in kwargs["measure"]:
        if l["genes"] == "all":
            x = range(X.shape[0])
        else:
            x = np.isin(adata.var[l["var"]].values, l["genes"])

        if "counts" == l:
            d["obs"][l["name"]] = np.array(X[:,x].sum(axis=1)).reshape(-1)
        elif "n_expressed_genes" == l:
            d["obs"][l["name"]] = np.array((X[:,x] > 0).sum(axis=1)).reshape(-1)
        elif "proportion":
            p = np.array(config.adata.X[:,x].sum(axis=1)).reshape(-1)/np.array(config.adata.X.sum(axis=1)).reshape(-1)
            d["obs"][l["name"]] =  p

    return d

def plot_violin(x=None, y=None, bandwidth=1):

    return go.Violin(
                        x=x,
                        y=y,
                        # bandwidth=bandwidth,
            )

def make_Graph(data=[],title="",x_label="",y_label="",style={"width": "90vw", "height": "50vh"}, center=True):
     
    graph =  dcc.Graph(
        figure={
                    "data":data,
                    "layout":{
                            "title": title,
                            "xaxis": {"title": x_label},
                            "yaxis": {"title": y_label},
                    }
        },
        style=style
    )

    if center:

        graph = [
            dbc.Col(),
            dbc.Col(graph),
            dbc.Col()
        ]

    return graph

def plot_qc():

    params = config.adata.uns[config.selected]["parameters"]

    if params["plot_style"] == "scatter":
        l = []
        # for i in data:

        #     x = None
        #     y = params["obs"][f"{config.selected}--{params['uns']['parameters']['plot_measure']}"].values

        #     l.append(
        #         make_Graph([plot_violin(x,y)], x_label=params['uns']['parameters']['batch'], y_label=params['uns']['parameters']['plot_measure'])
        #     )

        return l
        
    else:

        return []

    # node = get_node(name_analysis)
    # pos = get_node_pos(name_analysis)
    # if not node["data"]["computed"]:
    #     return []

    # dic = get_node_parameters(name_analysis,str2list=True)["measure"] 

    # node = get_node(name_analysis)
    # columns = node["data"]["threshold"]["columns"]
    # data = node["data"]["threshold"]["data"]

    # l = [
    #         html.H1("Thresholds"),
    #         dbc.Row(
    #             dash_table.DataTable(
    #                 id="qc_threshold_table",
    #                 columns=columns,
    #                 data=data
    #             )
    #         ),
    #     ]

    args = get_args(config.selected, adata)

    for metric in dic:
        var_selected_data = name_analysis+"_"+metric
            
        cols = [1]
        lims_min = get_table_column(data,metric+" min")
        lims_max = get_table_column(data,metric+" max")
        res = get_table_column(data,metric+" resolution")[0]

        custom_marker = {
            "symbol": "line-ns",  # Symbol code for a horizontal line (short dash)
            "size": 3,  # Length of the horizontal line
            "color": "red",  # Color of the line
            "line": {"width": 20}  # Line width
        }

        l += [
                dbc.Row([
                    dbc.Col(),
                    dbc.Col(
                        [dcc.Graph(id="Histogram",
                                figure={
                                        "data":[
                                            go.Scatter(
                                                x=cols,
                                                y=lims_min,
                                                mode="markers",
                                                marker=custom_marker,
                                                name="Local min threshold",
                                                # nbinsx=100,
                                                # name="Histogram",
                                                # opacity=0.7
                                            ),
                                            go.Scatter(
                                                x=cols,
                                                y=lims_max,
                                                mode="markers",
                                                marker=custom_marker,
                                                name="Local max threshold",
                                                # nbinsx=100,
                                                # name="Histogram",
                                                # marker=dict(color="blue"),
                                                # opacity=0.7
                                            ),
                                        ],
                                        "layout":{
                                                "title": f"Violinplot of {var_selected_data}",
                                                "xaxis": {"title": ""},
                                                "yaxis": {"title": var_selected_data},
                                                "barmode": "overlay",
                                        }
                                    },
                                style={"width": "90vw", "height": "50vh"}
                        )],
                    ),
                    dbc.Col(),
                    ],
                    justify="center"
                )
        ]
            
    lims_min_x = get_table_column(data,node["data"]["plot"]["x_label"]+" min")[0]
    lims_max_x = get_table_column(data,node["data"]["plot"]["x_label"]+" max")[0]
    lims_min_y = get_table_column(data,node["data"]["plot"]["y_label"]+" min")[0]
    lims_max_y = get_table_column(data,node["data"]["plot"]["y_label"]+" max")[0]

    l += [
            html.H1("Summary"),
            dcc.Dropdown(
                id="qc_summary_plot_x",
                value=node["data"]["plot"]["x_label"], 
                options=dic
            ),
            dcc.Dropdown(
                id="qc_summary_plot_y",
                value=node["data"]["plot"]["y_label"], 
                options=dic
            ),
            dcc.Dropdown(
                id="qc_summary_plot_c",
                value=node["data"]["plot"]["c_label"], 
                options=dic
            ),
            dbc.Row([
                dbc.Col(),
                dbc.Col(
                    [dcc.Graph(
                            id="summary_qc_scatter",
                            figure={
                                    "data":[
                                        go.Scatter(
                                            x=node["data"]["obs"][node["data"]["plot"]["x_label"]],
                                            y=node["data"]["obs"][node["data"]["plot"]["y_label"]],
                                            mode="markers",
                                            name="Local max threshold",
                                            # nbinsx=100,
                                            # name="Histogram",
                                            marker=dict(color=qualitative_colors(node["data"]["obs"][node["data"]["plot"]["c_label"]])),
                                            # opacity=0.7
                                        ),
                                        go.Scatter(
                                            x=[lims_min_x, lims_min_x, lims_max_x, lims_max_x, lims_min_x], 
                                            y=[lims_min_y,lims_max_y,lims_max_y,lims_min_y,lims_min_y],
                                            name="Limits",
                                            marker={"color":"black"}
                                        ),
                                    ],
                                    "layout":{
                                            "title": f"Summary plot",
                                            "xaxis": {"title": name_analysis+"_"+node["data"]["plot"]["x_label"]},
                                            "yaxis": {"title": name_analysis+"_"+node["data"]["plot"]["y_label"]},
                                            "barmode": "overlay",
                                    }
                                },
                            style={"width": "90vh", "height": "90vh"}
                    )],
                ),
                dbc.Col(),
                ],
                justify="center"
            )
    ]

    return l

# @app.callback(
#     dash.Output("analysis_plot","children",allow_duplicate=True),
#     dash.Input("qc_threshold_table","data"),
#     prevent_initial_call=True
# )
# def update_table(data):

#     prevent_race("qc")

#     node = get_node(config.selected)
#     pos = get_node_pos(config.selected)

#     dic = get_node_parameters(config.selected,str2list=True)["measure"] 
#     for i in dic:
#         col = np.array([int(i) for i in get_table_column(data,i+" resolution")])

#         s = np.array(node["data"]["obs"][i]) >= float(get_table_value(data,None,i+" min"))
#         s *= np.array(node["data"]["obs"][i]) <= float(get_table_value(data,None,i+" max"))
#         s = s < 1
#         config.graph[pos]["data"]["filter"][i] = s

#         v1 = mode(col).mode
#         if sum(col!=v1)>0:
#             v1 = col[col!=v1][0]

#         set_table_column(data,i+" resolution",v1)

#     pos = get_node_pos(config.selected)
#     config.graph[pos]["data"]["threshold"]["data"] = data

#     return plot_qc(config.selected)

# @app.callback(
#     dash.Output("summary_qc_scatter","figure",allow_duplicate=True),
#     dash.Input("qc_summary_plot_x","value"),
#     dash.Input("qc_summary_plot_y","value"),
#     dash.Input("qc_summary_plot_c","value"),
#     prevent_initial_call=True
# )
# def update_qc_summary(a,b,c):

#     prevent_race("qc")

#     node = get_node(config.selected)

#     pos = get_node_pos(config.selected)
#     node = get_node(config.selected)
#     data = node["data"]["thresholds"]["data"]
#     config.graph[pos]["data"]["x_label"] = a
#     config.graph[pos]["data"]["y_label"] = b
#     config.graph[pos]["data"]["c_label"] = c

#     lims_min_x = get_table_column(data,node["data"]["x_label"]+" min")[0]
#     lims_max_x = get_table_column(data,node["data"]["x_label"]+" max")[0]
#     lims_min_y = get_table_column(data,node["data"]["y_label"]+" min")[0]
#     lims_max_y = get_table_column(data,node["data"]["y_label"]+" max")[0]

#     f = {
#             "data":[
#                 go.Scatter(
#                     x=node["data"]["obs"][node["data"]["x_label"]],
#                     y=node["data"]["obs"][node["data"]["y_label"]],
#                     mode="markers",
#                     name="Local max threshold",
#                     # nbinsx=100,
#                     # name="Histogram",
#                     marker=dict(color=qualitative_colors(node["data"]["obs"][node["data"]["c_label"]])),
#                     # opacity=0.7
#                 ),
#                 go.Scatter(
#                     x=[lims_min_x, lims_min_x, lims_max_x, lims_max_x, lims_min_x], 
#                     y=[lims_min_y,lims_max_y,lims_max_y,lims_min_y,lims_min_y],
#                     name="Limits",
#                     marker={"color":"black"}
#                 ),
#             ],
#             "layout":{
#                     "title": f"Summary plot",
#                     "xaxis": {"title": config.selected+"_"+node["data"]["x_label"]},
#                     "yaxis": {"title": config.selected+"_"+node["data"]["y_label"]},
#                     "barmode": "overlay",
#             }
#         }

#     return f

config.methods["qc"] = {

    "properties":{"method":"qc","type":"QC","recompute":False},

    "args": args,

    "function":f_qc,

    "plot":plot_qc,
}