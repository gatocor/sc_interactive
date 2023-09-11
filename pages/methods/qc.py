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

    options = node_names(exclude_downstream_from_node=config.selected) 
    optionsthreshold = ["total_counts","n_genes_by_counts"]+[i for i in config.adata.var.columns.values if config.adata.var.dtypes[i] in [bool]]
    options_batch = get_batch_keys()

    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options":options
        },
        {
            "input":"Dropdown",
            "name":"batch_key",
            "description":"str, optional (default: None) Batch key to use. The Doublet metric will be computed independently in each set of cells separated by batch. If None, use the full dataset.",
            "value":None,
            "clearable":True,
            "options":options_batch,
            "summary":True
        },
        {
            "input":"MeasureTable",
            "name":"measure",
            "description":"Metrics of quality control to compute",
            "value":"[]",
            "options":optionsthreshold,
            "summary":True
        },
    ]

def f_qc(name_analysis, kwargs):
        
    for l in kwargs["measure"]:
        if "total_counts" == l:
            config.adata.obs[name_analysis+"_total_counts"] = np.array(config.adata.X.sum(axis=1)).reshape(-1)
        elif "n_genes_by_counts" == l:
            config.adata.obs[name_analysis+"_n_genes_by_counts"] = np.array((config.adata.X > 0).sum(axis=1)).reshape(-1)
        else:
            config.adata.obs[name_analysis+"_"+l] = np.array(config.adata.X[:,config.adata.var[l].values].sum(axis=1)).reshape(-1)/config.adata.obs["total_counts"].values
            config.adata.obs[name_analysis+"_"+l] = np.nan_to_num(config.adata.obs[l].values, nan=0)

    #Make empty table
    add = ["min","max","resolution"]
    if kwargs["batch_key"] == None:
        rows = [None]
    else:
        rows = np.sort(config.adata.obs[kwargs["batch_key"]].unique())
    
    columns, data = make_thresholds_table(kwargs["measure"], rows, add)
 
    #Fill table 
    for j in kwargs["measure"]:
        max_threshold = config.adata.obs[name_analysis+"_"+j].max()
        for i in rows:
            set_table_value(data, i, j+" min", 0)
            set_table_value(data, i, j+" max", max_threshold)
            set_table_value(data, i, j+" resolution", 100)

    pos = get_node_pos(config.selected)
    config.graph[pos]['data']['threshold'] = {'columns':columns,'data':data}

    if len(kwargs["measure"]) == 0:
        config.graph[pos]['data']['x_label'] = None
        config.graph[pos]['data']['y_label'] = None
        config.graph[pos]['data']['c_label'] = None
    elif len(kwargs["measure"]) == 1:
        config.graph[pos]['data']['x_label'] = kwargs["measure"][0]
        config.graph[pos]['data']['y_label'] = kwargs["measure"][0]
        config.graph[pos]['data']['c_label'] = kwargs["measure"][0]
    elif len(kwargs["measure"]) == 2:
        config.graph[pos]['data']['x_label'] = kwargs["measure"][0]
        config.graph[pos]['data']['y_label'] = kwargs["measure"][1]
        config.graph[pos]['data']['c_label'] = kwargs["measure"][0]
    elif len(kwargs["measure"]) >= 3:
        config.graph[pos]['data']['x_label'] = kwargs["measure"][0]
        config.graph[pos]['data']['y_label'] = kwargs["measure"][1]
        config.graph[pos]['data']['c_label'] = kwargs["measure"][2]

def rm_qc(name_analysis):

    ls = get_node_parameters(config.selected, str2list=True)["measure"]
    for l in ls:
        name = (name_analysis+"_"+l)
        config.adata.obs.drop(name, axis=1, inplace=True)

def rename_qc(name_analysis, name_new_analysis):

    ls = [(name_analysis+"_"+l) for l in get_node_parameters(name_analysis, str2list=True)["measure"]]
    cols = {}    
    for i in config.adata.obs.columns.values:
        if i not in ls:
            cols[i] = i
        else:
            name = i.split(name_analysis)[-1]
            name = (name_new_analysis+"_"+name)
            cols[i] = name

    config.adata.obs.rename(columns=cols, inplace=True)

def plot_qc(name_analysis):

    node = get_node(name_analysis)
    if not node['data']['computed']:
        return []

    dic = get_node_parameters(name_analysis,str2list=True)['measure'] 
    batch_key = get_node_parameters(name_analysis,str2list=True)['batch_key']

    node = get_node(name_analysis)
    columns = node['data']['threshold']['columns']
    data = node['data']['threshold']['data']

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
            
        if batch_key == None:

            cols = get_table_column(data,"Batch")

            lims_min = float(get_table_column(data,metric+" min")[0])
            lims_max = float(get_table_column(data,metric+" max")[0])
            res = int(get_table_column(data,metric+" resolution")[0])

            plot_max = hist_vline(config.adata.obs[var_selected_data].values, res)

            l += [
                    dbc.Row([
                        dbc.Col(),
                        dbc.Col(
                            [
                                dcc.Graph(id="Histogram",
                                    figure={
                                            "data":[
                                                go.Histogram(
                                                    x=config.adata.obs[var_selected_data].values,
                                                    nbinsx=res,
                                                    name='Histogram',
                                                    marker=dict(color='blue'),
                                                    opacity=0.7
                                                ),
                                                go.Scatter(
                                                    x=[lims_min, lims_min],
                                                    y=[0,plot_max],
                                                    name='Min threshold',
                                                    marker=dict(color='orange'),
                                                    opacity=0.7
                                                ),
                                                go.Scatter(
                                                    x=[lims_max, lims_max],
                                                    y=[0,plot_max],
                                                    name='Max threshold',
                                                    marker=dict(color='red'),
                                                    opacity=0.7
                                                ),
                                            ],
                                            "layout":{
                                                    'title': f'Histogram of {var_selected_data}',
                                                    'xaxis': {'title': var_selected_data},
                                                    'yaxis': {'title': 'Count'},
                                                    'barmode': 'overlay',
                                            }
                                        },
                                    style={'width': '90vw', 'height': '50vh'}
                                )
                            ],
                        ),
                        dbc.Col(),
                        ],
                        justify='center'
                    )
            ]
        else:

            cols = get_table_column(data,"Batch")
            lims_min = get_table_column(data,metric+" min")
            lims_max = get_table_column(data,metric+" max")
            res = get_table_column(data,metric+" resolution")[0]

            custom_marker = {
                'symbol': 'line-ns',  # Symbol code for a horizontal line (short dash)
                'size': 3,  # Length of the horizontal line
                'color': 'red',  # Color of the line
                'line': {'width': 20}  # Line width
            }

            x = config.adata.obs[batch_key].values

            l += [
                    dbc.Row([
                        dbc.Col(),
                        dbc.Col(
                            [dcc.Graph(id="Histogram",
                                    figure={
                                            "data":[
                                                go.Violin(
                                                    x=x,
                                                    y=config.adata.obs[var_selected_data].values,
                                                    bandwidth=res,
                                                    name='var_selected_data',
                                                    # marker=dict(color='blue'),
                                                    # opacity=0.7
                                                ),
                                                go.Scatter(
                                                    x=cols,
                                                    y=lims_min,
                                                    mode='markers',
                                                    marker=custom_marker,
                                                    name='Local min threshold',
                                                    # nbinsx=100,
                                                    # name='Histogram',
                                                    # opacity=0.7
                                                ),
                                                go.Scatter(
                                                    x=cols,
                                                    y=lims_max,
                                                    mode='markers',
                                                    marker=custom_marker,
                                                    name='Local max threshold',
                                                    # nbinsx=100,
                                                    # name='Histogram',
                                                    # marker=dict(color='blue'),
                                                    # opacity=0.7
                                                ),
                                            ],
                                            "layout":{
                                                    'title': f'Violinplot of {var_selected_data}',
                                                    'xaxis': {'title': batch_key},
                                                    'yaxis': {'title': var_selected_data},
                                                    'barmode': 'overlay',
                                            }
                                        },
                                    style={'width': '90vw', 'height': '50vh'}
                            )],
                        ),
                        dbc.Col(),
                        ],
                        justify='center'
                    )
            ]
            
    if batch_key == None:

        lims_min_x = get_table_column(data,node['data']['x_label']+" min")[0]
        lims_max_x = get_table_column(data,node['data']['x_label']+" max")[0]
        lims_min_y = get_table_column(data,node['data']['y_label']+" min")[0]
        lims_max_y = get_table_column(data,node['data']['y_label']+" max")[0]

        l += [
                html.H1("Summary"),
                dcc.Dropdown(
                    id="qc_summary_plot_x",
                    value=node['data']['x_label'], 
                    options=dic
                ),
                dcc.Dropdown(
                    id="qc_summary_plot_y",
                    value=node['data']['y_label'], 
                    options=dic
                ),
                dcc.Dropdown(
                    id="qc_summary_plot_c",
                    value=node['data']['c_label'], 
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
                                                x=config.adata.obs[name_analysis+"_"+node['data']['x_label']].values,
                                                y=config.adata.obs[name_analysis+"_"+node['data']['y_label']].values,
                                                mode='markers',
                                                name='Local max threshold',
                                                # nbinsx=100,
                                                # name='Histogram',
                                                marker=dict(color=qualitative_colors(config.adata.obs[name_analysis+"_"+node['data']['c_label']].values)),
                                                # opacity=0.7
                                            ),
                                            go.Scatter(
                                                x=[lims_min_x, lims_min_x, lims_max_x, lims_max_x, lims_min_x], 
                                                y=[lims_min_y,lims_max_y,lims_max_y,lims_min_y,lims_min_y],
                                                name='Limits',
                                                marker={'color':'black'}
                                            ),
                                        ],
                                        "layout":{
                                                'title': f'Summary plot',
                                                'xaxis': {'title': name_analysis+"_"+node['data']['x_label']},
                                                'yaxis': {'title': name_analysis+"_"+node['data']['y_label']},
                                                'barmode': 'overlay',
                                        }
                                    },
                                style={'width': '90vh', 'height': '90vh'}
                        )],
                    ),
                    dbc.Col(),
                    ],
                    justify='center'
                )
        ]

    return l

@app.callback(
    dash.Output("analysis_plot","children",allow_duplicate=True),
    dash.Input("qc_threshold_table","data"),
    prevent_initial_call=True
)
def update_table(data):

    prevent_race('qc')

    node = get_node(config.selected)

    dic = get_node_parameters(config.selected,str2list=True)['measure'] 
    batch = get_node_parameters(config.selected,str2list=True)['batch_key'] 
    for i in dic:
        if batch == None:
            col = np.array([int(i) for i in get_table_column(data,i+" resolution")])
        elif batch != None:
            col = np.array([float(i) for i in get_table_column(data,i+" resolution")])

        v1 = mode(col).mode
        if sum(col!=v1)>0:
            v1 = col[col!=v1][0]

        set_table_column(data,i+" resolution",v1)

    pos = get_node_pos(config.selected)
    config.graph[pos]['data']['threshold']['data'] = data

    return plot_qc(config.selected)

@app.callback(
    dash.Output("summary_qc_scatter","figure",allow_duplicate=True),
    dash.Input("qc_summary_plot_x","value"),
    dash.Input("qc_summary_plot_y","value"),
    dash.Input("qc_summary_plot_c","value"),
    prevent_initial_call=True
)
def update_qc_summary(a,b,c):

    prevent_race('qc')

    node = get_node(config.selected)

    pos = get_node_pos(config.selected)
    node = get_node(config.selected)
    data = node['data']['thresholds']['data']
    config.graph[pos]['data']['x_label'] = a
    config.graph[pos]['data']['y_label'] = b
    config.graph[pos]['data']['c_label'] = c

    lims_min_x = get_table_column(data,node['data']['x_label']+" min")[0]
    lims_max_x = get_table_column(data,node['data']['x_label']+" max")[0]
    lims_min_y = get_table_column(data,node['data']['y_label']+" min")[0]
    lims_max_y = get_table_column(data,node['data']['y_label']+" max")[0]

    f = {
            "data":[
                go.Scatter(
                    x=config.adata.obs[config.selected+"_"+node['data']['x_label']].values,
                    y=config.adata.obs[config.selected+"_"+node['data']['y_label']].values,
                    mode='markers',
                    name='Local max threshold',
                    # nbinsx=100,
                    # name='Histogram',
                    marker=dict(color=qualitative_colors(config.adata.obs[config.selected+"_"+node['data']['c_label']].values)),
                    # opacity=0.7
                ),
                go.Scatter(
                    x=[lims_min_x, lims_min_x, lims_max_x, lims_max_x, lims_min_x], 
                    y=[lims_min_y,lims_max_y,lims_max_y,lims_min_y,lims_min_y],
                    name='Limits',
                    marker={'color':'black'}
                ),
            ],
            "layout":{
                    'title': f'Summary plot',
                    'xaxis': {'title': config.selected+"_"+node['data']['x_label']},
                    'yaxis': {'title': config.selected+"_"+node['data']['y_label']},
                    'barmode': 'overlay',
            }
        }

    return f