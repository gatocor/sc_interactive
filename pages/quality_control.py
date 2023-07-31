import os
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import scanpy as sc
import dash_table
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from . import config
from .functions import *
import plotly.express as px
import scrublet

from app import app

click_doublet_old = 0

#Lists
def layout():
    return dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(html.H1("Quality control metrics"), width="auto"),
            ],
            justify="center",
            className="mb-4",
        ),
        dbc.Row(
                dash_table.DataTable(
                    id='table_qc_metrics',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Name"]
                    ],
                    data=[{"Name":i} for i in config.adata.uns["qc"]],
                    editable=False,
                    row_deletable=True,
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                ),
                # style={"margin-bottom":"1cm"}        
            ),
        dbc.Row(
            [
                dbc.Col(dcc.Dropdown(
                    id = "dropdown_add_metric",
                    value = None,
                    options = [str(i) for i in config.adata.var.columns.values if (config.adata.var.dtypes[i] in [bool]) and (i not in config.adata.uns["qc"])]
                )),
                dbc.Col(dbc.Button("Add metric", id="add_metric")),
            ],
            justify="Left",
            className="mb-4",
        ),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(html.H1("Global Quality control"), width="auto"),
            ],
            justify="center",
            className="mb-4",
        ),
        dbc.Row(
            [
                dbc.Col(html.H1("Thresholds"), width="auto"),
            ],
            justify="left",
            className="mb-4",
        ),
        dbc.Row(
            dash_table.DataTable(
                    id='qc_thresholds_table',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": True if ("Min" in i) or ("Max" in i) else False} for i in config.qc_df_threshold.columns
                    ],
                    data=f_qc_table_threshold(config.adata),
                    editable=False,
                    row_deletable=False,
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                ),
                style={"margin-bottom":"1cm"}
        ),
        dbc.Row(
            [
                dbc.Col(html.H1("Histogram plots"), width="auto"),
            ],
            justify="left",
            className="mb-4",
        ),
        dbc.Row(
            id = "global_qc_plots",
            children = plot_qc_global_hist(config.adata),
            justify="center",
            className="mb-4"
        ),
        dbc.Row(
            [
                dbc.Col(html.H1("Summary plot"), width="auto"),
            ],
            justify="left",
            className="mb-4",
        ),
        dbc.Row(
            [
                dbc.Row([
                    dbc.Col(
                        [
                            dbc.Tooltip(
                                "X axis",
                                target="y-size-label",
                                placement="bottom",
                            ),
                            html.Label("Plot style:", className="input-label", id="y-size-label"),
                        ],
                        width=2,
                    ),
                    dbc.Col(
                        [
                            html.Div(
                                dcc.Dropdown(
                                    id='qc_scatter1_dropdown',
                                    options=[i for i in config.adata.uns["qc"]],
                                    value=config.qc_summary_x,
                                    placeholder="Select a column",
                                    clearable=False
                                )
                            ),
                        ]
                    ),
                    ]
                ),
                dbc.Row([
                    dbc.Col(
                        [
                            dbc.Tooltip(
                                "Y axis",
                                target="y-size-label",
                                placement="bottom",
                            ),
                            html.Label("Plot style:", className="input-label", id="y-size-label"),
                        ],
                        width=2,
                    ),
                    dbc.Col(
                        [
                            html.Div(
                                dcc.Dropdown(
                                    id='qc_scatter2_dropdown',
                                    options=[i for i in config.adata.uns["qc"]],
                                    value=config.qc_summary_y,
                                    placeholder="Select a column",
                                    clearable=False
                                )
                            ),
                        ]
                    ),
                    ]
                ),
                dbc.Row([
                    dbc.Col(
                        [
                            dbc.Tooltip(
                                "Color",
                                target="y-size-label",
                                placement="bottom",
                            ),
                            html.Label("Plot style:", className="input-label", id="y-size-label"),
                        ],
                        width=2,
                    ),
                    dbc.Col(
                        [
                            html.Div(
                                dcc.Dropdown(
                                    id='qc_scatter3_dropdown',
                                    options=[i for i in config.adata.obs.columns.values],
                                    value=config.qc_summary_color,
                                    placeholder="Select a column",
                                    clearable=False
                                )
                            ),
                        ]
                    ),
                    ]
                ),
                dbc.Row(
                    dbc.Col(
                            dcc.Graph(id='qc_plot_scatter')
                    ),
                    style={'width': '50%', 'margin': 'auto'}
                )
            ]
        ),
        html.Hr(),
        dbc.Row(
            id = 'qc_per_condition',
            children = make_qc_per_condition(config.adata)
            ,
            justify="center",
            className="mb-4"
        ),
        html.Hr(),
        dbc.Row(
            id = 'doublets',
            children = make_qc_doublets(config.adata),
            justify="center",
            className="mb-4"
        ),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(html.H1("Summary"), width="auto"),
            ],
            justify="center",
            className="mb-4"
        ),
        dbc.Row(
                dash_table.DataTable(
                    id='qc_summary_table',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Concept","Removed cells","Removed cells (%)"]
                    ],
                    data=qc_summary(config.adata),
                    editable=False,
                    row_deletable=True,
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                ),
                style={"margin-bottom":"1cm"}        
            ),
        dbc.Row(
            [
                dbc.Button(id='save_qc-button', n_clicks=0, children="Save",
                            size="lg",
                            style={
                                "background-color":"#343A40",
                                   'width': '280px', 
                            }      
                            )
            ],
            justify="center",
            className="mb-4"
        )
    ],
    fluid=True
)

@app.callback(
    dash.Output('table_qc_metrics','data'),
    dash.Output('qc_thresholds_table', 'data'),
    dash.Output("global_qc_plots","children"),
    dash.Output("dropdown_add_metric","options"),
    dash.Output("dropdown_add_metric","value"),
    dash.Output('qc_scatter1_dropdown', 'options'),
    dash.Output('qc_scatter2_dropdown', 'options'),
    dash.Output('qc_scatter3_dropdown', 'options'),
    dash.Output('qc_summary_table', 'data', allow_duplicate=True),
    [
        dash.Input('add_metric', 'n_clicks'),
        dash.Input('table_qc_metrics', 'data'),
        dash.Input('qc_thresholds_table', 'data'),
    ],
    [
        dash.State('dropdown_add_metric', 'value'),
    ],
    prevent_initial_call=True
)
def update_qc_global_histogram(n_clicks, qc_metrics, data, value):
    
    if value != None:
        qc_metrics.append({"Name":str(value)})

    f_qc(config.adata, qc_metrics, data)
        
    options = [str(i) for i in config.adata.var.columns.values if (config.adata.var.dtypes[i] in [bool]) and (i not in config.adata.uns["qc"])]
        
    options_dropdown = list(config.adata.uns["qc"].keys())
    options_dropdown_color = [i for i in config.adata.obs.columns.values]

    return f_qc_table_metrics(config.adata), f_qc_table_threshold(config.adata), plot_qc_global_hist(config.adata), options, None, options_dropdown, options_dropdown, options_dropdown_color, qc_summary(config.adata) #{'data': data,'layout': layout}

@app.callback(
    dash.dependencies.Output('qc_plot_scatter', 'figure'),
    [
     dash.Input('qc_scatter1_dropdown', 'value'),
     dash.Input('qc_scatter2_dropdown', 'value'),
     dash.Input('qc_scatter3_dropdown', 'value'),
     dash.Input('qc_thresholds_table', 'data')        
     ],
)
def update_qc_global_scatter(var_selected_data1, var_selected_data2, var_selected_data3, data):
    
    config.qc_summary_x = var_selected_data1
    config.qc_summary_y = var_selected_data2
    config.qc_summary_color = var_selected_data3

    # Create a vertical line at the specified input value
    var1_vertical_line_min = go.Scatter(
        x=[config.adata.uns["qc"][var_selected_data1]["Minimum threshold"], config.adata.uns["qc"][var_selected_data1]["Minimum threshold"]],
        y=[np.min(config.adata.obs[var_selected_data2].values), np.max(config.adata.obs[var_selected_data2].values)],
        mode='lines',
        name='Min threshold',
        line=dict(color='red', dash='dash')
    )
    var1_vertical_line_max = go.Scatter(
        x=[config.adata.uns["qc"][var_selected_data1]["Maximum threshold"], config.adata.uns["qc"][var_selected_data1]["Maximum threshold"]],
        y=[np.min(config.adata.obs[var_selected_data2].values), np.max(config.adata.obs[var_selected_data2].values)],
        mode='lines',
        name='Max threshold',
        line=dict(color='green', dash='dash')
    )

    # Create a vertical line at the specified input value
    var2_vertical_line_min = go.Scatter(
        y=[config.adata.uns["qc"][var_selected_data2]["Minimum threshold"], config.adata.uns["qc"][var_selected_data2]["Minimum threshold"]],
        x=[np.min(config.adata.obs[var_selected_data1].values), np.max(config.adata.obs[var_selected_data1].values)],
        mode='lines',
        name='Min threshold',
        line=dict(color='red', dash='dot')
    )
    var2_vertical_line_max = go.Scatter(
        y=[config.adata.uns["qc"][var_selected_data2]["Maximum threshold"], config.adata.uns["qc"][var_selected_data2]["Maximum threshold"]],
        x=[np.min(config.adata.obs[var_selected_data1].values), np.max(config.adata.obs[var_selected_data1].values)],
        mode='lines',
        name='Max threshold',
        line=dict(color='green', dash='dot')
    )

    if config.adata.obs.dtypes[var_selected_data3] in ["category", object, int, str]:
        labels = config.adata.obs[var_selected_data3].values
        # Use the default Plotly Express color sequence
        color_sequence = px.colors.qualitative.Plotly

        # Create a dictionary to map categories to colors
        color_map = {
            category: color_sequence[i % len(color_sequence)]
            for i, category in enumerate(np.unique(labels))
        }
        sc = go.Scatter(
            x= config.adata.obs[var_selected_data1].values,
            y= config.adata.obs[var_selected_data2].values,
            marker=dict(
                color=[color_map[i] for i in labels],
            ),
            mode='markers',
            name='Scatter',
        )
    else:
        sc = go.Scatter(
            x= config.adata.obs[var_selected_data1].values,
            y= config.adata.obs[var_selected_data2].values,
            marker=dict(
                color=config.adata.obs[var_selected_data3].values,
            ),
            mode='markers',
            name='Scatter',
        )

    data = [
        sc,
        var1_vertical_line_min,
        var1_vertical_line_max,
        var2_vertical_line_min,
        var2_vertical_line_max,
    ]

    layout = {
        'xaxis': {'title': var_selected_data1},
        'yaxis': {'title': var_selected_data2,
                    # 'scaleanchor': 'x',
                    'scaleratio': 1
                  },
        'barmode': 'overlay',
        'width':900,
        'height':800,
    }

    return {
        'data': data,
        'layout': layout,
    }

@app.callback(
     dash.Output('qc_plot_violin', 'figure'),
    [
        dash.Input('qc_per_condition_dropdown2', 'value'),
        dash.Input('qc_per_condition_dropdown1', 'value'),
        dash.Input('qc_thresholds_table', 'data'),
    ]
)
def update_qc_threshold_table(var_selected_data1, var_selected_data2, data):

    #Plot_type
    # Create a vertical line at the specified input value
    var1_vertical_line_min = go.Scatter(
        x=config.adata.obs[var_selected_data2].values,
        y=[config.adata.uns["qc"][var_selected_data1]["Minimum threshold"] for i in range(len(config.adata.obs[var_selected_data2].values))],
        mode='lines',
        name='Min threshold',
        line=dict(color='red')
    )
    var1_vertical_line_max = go.Scatter(
        x=config.adata.obs[var_selected_data2].values,
        y=[config.adata.uns["qc"][var_selected_data1]["Maximum threshold"] for i in range(len(config.adata.obs[var_selected_data2].values))],
        mode='lines',
        name='Max threshold',
        line=dict(color='green')
    )

    order = np.argsort(config.adata.obs[var_selected_data2].values)

    data = [
        go.Violin(
            x=config.adata.obs[var_selected_data2].values[order],
            y=config.adata.obs[var_selected_data1].values[order],
            name='Violin',
            # marker=dict(color='blue'),
            opacity=0.7
        ),
        var1_vertical_line_min,
        var1_vertical_line_max
    ]

    layout = {
        'title': f'Violinplot of {var_selected_data1} per {var_selected_data2}',
        'xaxis': {'title': var_selected_data2},
        'yaxis': {'title': var_selected_data1},
        'barmode': 'overlay',
        # 'width':900,
        'height':800,
    }

    return {
        'data': data,
        'layout': layout,
    }

@app.callback(
     dash.Output('qc_thresholds_per_condition_table', 'columns'),
     dash.Output('qc_thresholds_per_condition_table', 'data'),
    [
        dash.Input('qc_per_condition_dropdown1', 'value'),
        dash.Input('qc_per_condition_dropdown2', 'value'),
    ]
)
def update_qc_threshold_per_condition_table(condition,qc):

    data = []

    if isinstance(config.adata.uns["qc"][qc]["Minimum threshold"],list):
        data.append(
            {j:config.adata.uns["qc"][qc]["Minimum threshold"][i] for i,j in enumerate(config.adata.obs[condition].cat.categories)}
        )

        data.append(
            {j:config.adata.uns["qc"][qc]["Maximum threshold"][i] for i,j in enumerate(config.adata.obs[condition].cat.categories)}
        )
    else:
        data.append(
            {j:config.adata.uns["qc"][qc]["Minimum threshold"] for i,j in enumerate(config.adata.obs[condition].cat.categories)}
        )

        data.append(
            {j:config.adata.uns["qc"][qc]["Maximum threshold"] for i,j in enumerate(config.adata.obs[condition].cat.categories)}
        )

    cols = [{"name": i, "id": i, "deletable": False, "editable": True} for i in config.adata.obs[condition].cat.categories]

    return cols, data 

@app.callback(
     dash.Output('qc_per_condition', 'children'),
    [
        dash.Input('qc_per_condition-button', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def make_qc_per_condition_(n_clicks):

    if "per_condition" in config.adata.uns["qc"]["total_counts"].keys():

        for var_selected_data in config.adata.uns["qc"].keys():          
            del config.adata.uns["qc"][var_selected_data]["per_condition"]

    else:

        for var_selected_data in config.adata.uns["qc"].keys():          
            config.adata.uns["qc"][var_selected_data]["per_condition"] = {}

    return make_qc_per_condition(config.adata)

@app.callback(
     dash.Output('per_condition_plot', 'children'),
     dash.Output('table_qc_per_condition_metrics', 'data'),
     dash.Output('dropdown_add_per_condition_metrics', 'value'),
     dash.Output('table_qc_per_condition', 'data'),
     dash.Output('qc_summary_table', 'data', allow_duplicate=True),
    [
        dash.Input('add_qc_per_condition-button', 'n_clicks'),
        dash.Input('table_qc_per_condition_metrics', 'data'),
        dash.Input('table_qc_per_condition', 'data'),
    ],
    dash.State('dropdown_add_per_condition_metrics', 'value'),
    prevent_initial_call=True
)
def update_qc_threshold_per_condition_table(condition, data, local_thresholds, value):

    conditions = [i["Name"] for i in data]
    if value != None and value not in conditions:
        data.append({"Name":str(value)})

    conditions = [i["Name"] for i in data]

    for var_selected_data in [i for i in config.adata.uns["qc"]]:          
        
        if len(data) != 0:                
            for condition in conditions:
                if condition not in config.adata.uns["qc"][var_selected_data]["per_condition"].keys():
                    config.adata.uns["qc"][var_selected_data]["per_condition"][condition] = \
                        {
                            str(ci):{
                                "Min":float(config.adata.uns["qc"][var_selected_data]["Minimum threshold"]),
                                "Max":float(config.adata.uns["qc"][var_selected_data]["Maximum threshold"])                                                              
                            }
                            for ci in np.unique(config.adata.obs[condition].values)
                        }

    for i in local_thresholds:
        if i["Min"] == None:
            config.adata.uns["qc"][i["Variable"]]["per_condition"][i["Condition"]][i["Condition type"]]["Min"] = config.adata.uns["qc"][i["Variable"]]["Minimum threshold"]
        else:
            config.adata.uns["qc"][i["Variable"]]["per_condition"][i["Condition"]][i["Condition type"]]["Min"] = float(i["Min"])

        if i["Max"] == None:
            config.adata.uns["qc"][i["Variable"]]["per_condition"][i["Condition"]][i["Condition type"]]["Max"] = config.adata.uns["qc"][i["Variable"]]["Maximum threshold"]
        else:
            config.adata.uns["qc"][i["Variable"]]["per_condition"][i["Condition"]][i["Condition type"]]["Max"] = float(i["Max"])
            
    return qc_per_condition_plots(config.adata), data, None, qc_per_condition_table(config.adata), qc_summary(config.adata)

@app.callback(
     dash.Output('doublets', 'children'),
    [
        dash.Input('doublets-button', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def make_qc_doublets_(n_clicks):

    if "doublets" in config.adata.uns.keys():

        del config.adata.uns["doublets"]
    
    else:

        config.adata.uns["doublets"] = {}

    return make_qc_doublets(config.adata)

@app.callback(
     dash.Output('qc_plot_doublets', 'children', allow_duplicate=True),
     dash.Output('table_doublets_threshold', 'data', allow_duplicate=True),
     dash.Output('qc_summary_table', 'data', allow_duplicate=True),
    [
        dash.Input('button_doublets', 'n_clicks'),
    ],
    [
        dash.State('doublets_batch_key', 'value'),
        dash.State('doublets_qc_before_computation', 'value'),
        dash.State('doublets_synthetic_doublet_umi_subsampling', 'value'),
        dash.State('doublets_use_approx_neighbors', 'value'),
        dash.State('doublets_distance_metric', 'value'),
        dash.State('doublets_min_counts', 'value'),
        dash.State('doublets_min_cells', 'value'),
        dash.State('doublets_min_gene_variability_pctl', 'value'),
        dash.State('doublets_log_transform', 'value'),
        dash.State('doublets_mean_center', 'value'),
        dash.State('doublets_normalize_variance', 'value'),
        dash.State('doublets_n_prin_comps', 'value'),
        dash.State('doublets_svd_solver', 'value'),
    ],
    prevent_initial_call=True
)
def compute_qc_doublets(
    n_clicks, 
    batch_key, 
    qc_before_computation, 
    synthetic_doublet_umi_subsampling, 
    use_approx_neighbors, 
    distance_metric, 
    min_counts, min_cells, 
    min_gene_variability_pctl, 
    log_transform, 
    mean_center, 
    normalize_variance, 
    n_prin_comps, 
    svd_solver):

    l = []
    if n_clicks != None:

        scrub = scrublet.Scrublet(config.adata.X)

        scrub.scrub_doublets(
                                synthetic_doublet_umi_subsampling = synthetic_doublet_umi_subsampling, 
                                use_approx_neighbors = use_approx_neighbors, 
                                distance_metric = distance_metric, 
                                min_counts = min_counts,
                                min_cells = min_cells, 
                                min_gene_variability_pctl = min_gene_variability_pctl, 
                                log_transform = log_transform, 
                                mean_center = mean_center, 
                                normalize_variance = normalize_variance, 
                                n_prin_comps = n_prin_comps, 
                                svd_solver= svd_solver
                            )
            
        X = scrublet.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)

        threshold = {"full":1}
        if batch_key == "full":
            threshold = {i:1 for i in np.unique(config.adata.obs[batch_key].values)}

        config.adata.obs["doublet_score"] = scrub.doublet_scores_obs_
        config.adata.uns["doublets"] = \
            {
                "algorithm":"scrublet",
                "doublet_score_sim": scrub.doublet_scores_sim_,
                "X_umap": X,
                "batch_key": batch_key,
                "threshold": threshold,
                "params":{
                    "synthetic_doublet_umi_subsampling": synthetic_doublet_umi_subsampling, 
                    "use_approx_neighbors": use_approx_neighbors, 
                    "distance_metric": distance_metric, 
                    "min_counts": min_counts,
                    "min_cells": min_cells, 
                    "min_gene_variability_pctl": min_gene_variability_pctl, 
                    "log_transform": log_transform, 
                    "mean_center": mean_center, 
                    "normalize_variance": normalize_variance, 
                    "n_prin_comps": n_prin_comps, 
                    "svd_solver": svd_solver
                }
            }
        
        l = [{"Condition":i, "Threshold":j} for i,j in config.adata.uns["doublets"]["threshold"].items()]

    return make_qc_doublet_plots(config.adata), l, qc_summary(config.adata)

@app.callback(
     dash.Output('qc_plot_doublets', 'children', allow_duplicate=True),
     dash.Output('table_doublets_threshold', 'data', allow_duplicate=True),
     dash.Output('qc_summary_table', 'data', allow_duplicate=True),
    [
        dash.Input('table_doublets_threshold', 'data'),
    ],
    prevent_initial_call=True
)
def compute_qc_doublets(
    data
    ):

    l = []
    if data != []:
        config.adata.uns["doublets"]["threshold"] = {i["Condition"]:float(i["Threshold"]) for i in data}

        l = [{"Condition":i, "Threshold":j} for i,j in config.adata.uns["doublets"]["threshold"].items()]

    return make_qc_doublet_plots(config.adata), l, qc_summary(config.adata)

@app.callback(
     dash.Output('save_qc-button', 'children'),
    [
     dash.Input('save_qc-button', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def save_gene_lists(n_clicks):

    config.adata.write(config.file_path)

    l = qc_limit(config.adata)

    config.adata = config.adata[l].copy()

    config.adata.write(config.file_path.split(".h5ad")[0]+"_qc.h5ad")

    return "Save"