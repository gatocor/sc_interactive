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

from ..functions import *

from app import app

from .. import config

def args_qc():

    optionsthreshold = ["total_counts","n_genes_by_counts"]+[i for i in config.adata.var.columns.values if config.adata.var.dtypes[i] in [bool]]

    return [
        # {
        #     "input":"MeasureTable",
        #     "name":"measure",
        #     "description":"Metrics of quality control to compute",
        #     "value":"[total_counts,n_genes_by_counts]",
        #     "options":optionsthreshold,
        #     "summary":True
        # },
        {
            "input":"GeneTable",
            "name":"measure",
            "description":"Metrics of quality control to compute",
            "value":"[total_counts,n_genes_by_counts]",
            "options":optionsthreshold,
            "summary":True
        },
    ]

def f_qc(name_analysis, kwargs):
        
    pos = get_node_pos(name_analysis)

    config.graph[pos]["data"]["obs"] = {}
    for l in kwargs["measure"]:
        if "total_counts" == l:
            config.graph[pos]["data"]["obs"][l] = np.array(config.adata.X.sum(axis=1)).reshape(-1)
        elif "n_genes_by_counts" == l:
            config.graph[pos]["data"]["obs"][l] = np.array((config.adata.X > 0).sum(axis=1)).reshape(-1)
        else:
            config.graph[pos]["data"]["obs"][l] = np.array(config.adata.X[:,config.adata.var[l].values].sum(axis=1)).reshape(-1)/np.array(config.adata.X.sum(axis=1)).reshape(-1)
            config.graph[pos]["data"]["obs"][l] = np.nan_to_num(config.adata.obs[l].values, nan=0)

    #Make empty table
    add = ["min","max","resolution"]
    rows = [None]

    columns, data = make_thresholds_table(kwargs["measure"], rows, add)
 
    #Fill table 
    for j in kwargs["measure"]:
        max_threshold = np.max(config.graph[pos]["data"]["obs"][j])
        for i in rows:
            set_table_value(data, i, j+" min", 0)
            set_table_value(data, i, j+" max", max_threshold)
            set_table_value(data, i, j+" resolution", 100)

    pos = get_node_pos(config.selected)
    config.graph[pos]["data"]["threshold"] = {"columns":columns,"data":data}
    config.graph[pos]["data"]["filter"] = {i:np.ones_like(config.adata.X.shape[0])>0 for i in get_node_parameters(name_analysis,str2list=True)["measure"]}

    if len(kwargs["measure"]) == 0:
        config.graph[pos]["data"]["plot"] = {
            "x_label": None,
            "y_label": None,
            "c_label": None
        }
    elif len(kwargs["measure"]) == 1:
        config.graph[pos]["data"]["plot"] = {
            "x_label": kwargs["measure"][0],
            "y_label": kwargs["measure"][0],
            "c_label": kwargs["measure"][0],
        }
    elif len(kwargs["measure"]) == 2:
        config.graph[pos]["data"]["plot"] = {
            "x_label": kwargs["measure"][0],
            "y_label": kwargs["measure"][1],
            "c_label": kwargs["measure"][0],
        }
    elif len(kwargs["measure"]) >= 3:
        config.graph[pos]["data"]["plot"] = {
            "x_label": kwargs["measure"][0],
            "y_label": kwargs["measure"][1],
            "c_label": kwargs["measure"][2],
        }

def rm_qc(name_analysis):

    return

def rename_qc(name_analysis, name_new_analysis):

    return

def plot_qc(name_analysis):

    node = get_node(name_analysis)
    pos = get_node_pos(name_analysis)
    if not node["data"]["computed"]:
        return []

    dic = get_node_parameters(name_analysis,str2list=True)["measure"] 

    node = get_node(name_analysis)
    columns = node["data"]["threshold"]["columns"]
    data = node["data"]["threshold"]["data"]

    l = [
            html.H1("Thresholds"),
            dbc.Row(
                dash_table.DataTable(
                    id="qc_threshold_table",
                    columns=columns,
                    data=data
                )
            ),
        ]

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
                                            go.Violin(
                                                y=node["data"]["obs"][metric],
                                                bandwidth=res,
                                                name="var_selected_data",
                                                # marker=dict(color="blue"),
                                                # opacity=0.7
                                            ),
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

@app.callback(
    dash.Output("analysis_plot","children",allow_duplicate=True),
    dash.Input("qc_threshold_table","data"),
    prevent_initial_call=True
)
def update_table(data):

    prevent_race("qc")

    node = get_node(config.selected)
    pos = get_node_pos(config.selected)

    dic = get_node_parameters(config.selected,str2list=True)["measure"] 
    for i in dic:
        col = np.array([int(i) for i in get_table_column(data,i+" resolution")])

        s = np.array(node["data"]["obs"][i]) >= float(get_table_value(data,None,i+" min"))
        s *= np.array(node["data"]["obs"][i]) <= float(get_table_value(data,None,i+" max"))
        s = s < 1
        config.graph[pos]["data"]["filter"][i] = s

        v1 = mode(col).mode
        if sum(col!=v1)>0:
            v1 = col[col!=v1][0]

        set_table_column(data,i+" resolution",v1)

    pos = get_node_pos(config.selected)
    config.graph[pos]["data"]["threshold"]["data"] = data

    return plot_qc(config.selected)

@app.callback(
    dash.Output("summary_qc_scatter","figure",allow_duplicate=True),
    dash.Input("qc_summary_plot_x","value"),
    dash.Input("qc_summary_plot_y","value"),
    dash.Input("qc_summary_plot_c","value"),
    prevent_initial_call=True
)
def update_qc_summary(a,b,c):

    prevent_race("qc")

    node = get_node(config.selected)

    pos = get_node_pos(config.selected)
    node = get_node(config.selected)
    data = node["data"]["thresholds"]["data"]
    config.graph[pos]["data"]["x_label"] = a
    config.graph[pos]["data"]["y_label"] = b
    config.graph[pos]["data"]["c_label"] = c

    lims_min_x = get_table_column(data,node["data"]["x_label"]+" min")[0]
    lims_max_x = get_table_column(data,node["data"]["x_label"]+" max")[0]
    lims_min_y = get_table_column(data,node["data"]["y_label"]+" min")[0]
    lims_max_y = get_table_column(data,node["data"]["y_label"]+" max")[0]

    f = {
            "data":[
                go.Scatter(
                    x=node["data"]["obs"][node["data"]["x_label"]],
                    y=node["data"]["obs"][node["data"]["y_label"]],
                    mode="markers",
                    name="Local max threshold",
                    # nbinsx=100,
                    # name="Histogram",
                    marker=dict(color=qualitative_colors(node["data"]["obs"][node["data"]["c_label"]])),
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
                    "xaxis": {"title": config.selected+"_"+node["data"]["x_label"]},
                    "yaxis": {"title": config.selected+"_"+node["data"]["y_label"]},
                    "barmode": "overlay",
            }
        }

    return f