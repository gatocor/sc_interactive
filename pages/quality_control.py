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
from .args.doublet_args import doublet_args

from app import app

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
                    data=[{"Name":"total_counts"},{"Name":"n_genes_by_counts"}],
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
                    value = ([[i for i in config.adata.uns["gene_lists"].keys()]]+[None])[0],
                    options = [i for i in config.adata.uns["gene_lists"].keys()]
                )),
                dbc.Col(dbc.Button("Add metric")),
            ],
            justify="Left",
            className="mb-4",
        ),
        dbc.Row(
            [
                dbc.Col(html.H1("Global Quality control"), width="auto"),
            ],
            justify="center",
            className="mb-4",
        ),
        dbc.Row(
            dash_table.DataTable(
                    id='qc_thresholds_table',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": True if ("Min" in i) or ("Max" in i) else False} for i in config.qc_df_threshold.columns
                    ],
                    data=config.qc_df_threshold.to_dict('records'),
                    editable=False,
                    row_deletable=False,
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                ),
                style={"margin-bottom":"1cm"}
        ),
        dbc.Row(
            [
                dbc.Col([
                    dbc.Row([
                        dbc.Col(
                            [
                                dbc.Tooltip(
                                    "Plot style",
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
                                        id='qc_hist_dropdown',
                                        options=f_options(config.adata,"qc_"),
                                        value="qc_total_counts",
                                        placeholder="Select a column",
                                        clearable=False
                                    )
                                ),
                            ]
                        ),
                    ]),
                ],
                width=6
                ),
                dbc.Col(
                    [
                        dbc.Row([
                            dbc.Col(
                                [
                                    dbc.Tooltip(
                                        "Plot style",
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
                                            options=f_options(config.adata,"qc_"),
                                            value="qc_total_counts",
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
                                        "Plot style",
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
                                            options=f_options(config.adata,"qc_"),
                                            value="qc_n_genes_by_counts",
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
                                        "Plot style",
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
                                            options=f_options(config.adata,"qc_"),
                                            value="qc_total_counts",
                                            placeholder="Select a column",
                                            clearable=False
                                        )
                                    ),
                                ]
                            ),
                            ]
                        ),
                    ],
                    width=6
                ),
            ],
            justify="center",
            align="top"
        ),
        dbc.Row([
            dbc.Col(
                        dcc.Graph(id='qc_plot_histogram'),
            ),
            dbc.Col(
                        dcc.Graph(id='qc_plot_scatter'),
            )
        ]),
        html.Hr(),
        dbc.Row(
            [
                dbc.Col(html.H1("Per Condition Quality control"), width="auto"),
            ],
            justify="center",
            className="mb-4"
        ),
        dbc.Row(
            dcc.Dropdown(
                id='qc_per_condition_dropdown1',
                options=[i for i in config.adata.obs.columns.values if i.startswith("condition_")],
                value=([i for i in config.adata.obs.columns.values if i.startswith("condition_")]+[0])[0],
                placeholder="Select a column",
                clearable=False
            )
        ),
        dbc.Row(
            dcc.Dropdown(
                id='qc_per_condition_dropdown2',
                options=[i for i in config.adata.obs.columns.values if i.startswith("qc_")],
                value=[i for i in config.adata.obs.columns.values if i.startswith("qc_")][0],
                placeholder="Select a column",
                clearable=False
            )
        ),
        dbc.Row(
                    dcc.Graph(id='qc_plot_violin'),
        ),  
        dash_table.DataTable(
                id='qc_thresholds_per_condition_table',
                columns=[
                    {"name": i, "id": i, "deletable": False, "editable": True if ("Min" in i) or ("Max" in i) else False} for i in config.qc_df_threshold.columns
                ],
                data=config.qc_df_threshold.to_dict('records'),
                editable=False,
                row_deletable=False,
                style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
            ),
        html.Hr(),
        dbc.Row(
            id = "data",
            children = [
                dbc.Button(id='doublets-button', n_clicks=0, children="Add Doublet Removal",
                            size="lg",
                            style={
                                "background-color":"#343A40",
                                   'width': '280px', 
                            }      
                            ),
            ],
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
            [
                dcc.Graph(
                    id='qc_summary_table',
                    figure=go.Figure(
                        data=[go.Table(
                            header=dict(values=config.qc_summary.columns),
                            cells=dict(values=[config.qc_summary[col] for col in config.qc_summary.columns])
                        )],
                    )
                )
            ],
            justify="center",
            className="mb-4"
        ),
        dbc.Row(
            [
                dbc.Button(id='save-button', n_clicks=0, children="Save",
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
    dash.Output('qc_plot_histogram', 'figure'),
    [
     dash.Input('qc_hist_dropdown', 'value'),
     dash.Input('qc_thresholds_table', 'data')        
     ],
)
def update_qc_global_histogram(var_selected_data, data):
    
    #Plot_type
    # Create a vertical line at the specified input value
    hist_values, hist_bins = np.histogram(config.adata.obs[var_selected_data].values, bins=30)
    tallest_bin_height = np.max(hist_values)
    var1_vertical_line_min = go.Scatter(
        x=[config.adata.uns["qc"][var_selected_data]["Minimum threshold"], config.adata.uns["qc"][var_selected_data]["Minimum threshold"]],
        y=[0, tallest_bin_height],
        mode='lines',
        name='Min threshold',
        line=dict(color='red', dash='dash')
    )
    var1_vertical_line_max = go.Scatter(
        x=[config.adata.uns["qc"][var_selected_data]["Maximum threshold"], config.adata.uns["qc"][var_selected_data]["Maximum threshold"]],
        y=[0, tallest_bin_height],
        mode='lines',
        name='Max threshold',
        line=dict(color='green', dash='dash')
    )

    data = [
        go.Histogram(
            x=config.adata.obs[var_selected_data].values,
            nbinsx=30,
            name='Histogram',
            marker=dict(color='blue'),
            opacity=0.7
        ),
        var1_vertical_line_min,
        var1_vertical_line_max
    ]

    layout = {
        'title': f'Histogram of {var_selected_data}',
        'xaxis': {'title': var_selected_data},
        'yaxis': {'title': 'Count'},
        'barmode': 'overlay',
        'width':900,
        'height':800,
    }


    return {
        'data': data,
        'layout': layout,
    }

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

    data = [
        go.Scatter(
            x= config.adata.obs[var_selected_data1].values,
            y= config.adata.obs[var_selected_data2].values,
            marker=dict(
                color=config.adata.obs[var_selected_data3].values,
            ),
            mode='markers',
            name='Scatter',
        ),
        var1_vertical_line_min,
        var1_vertical_line_max,
        var2_vertical_line_min,
        var2_vertical_line_max,
    ]

    layout = {
        'title': f'Scatter of {var_selected_data1} vs {var_selected_data2}',
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

# # Step 5: Run the app
# if __name__ == '__main__':
#     app.run_server(debug=True)

# Callback to handle adding and deleting rows
# @app.callback(
#      dash.Output('qc_pattern_table', 'data'),
#      dash.Output('qc_thresholds_table', 'data'),
#      dash.Output('qc_hist_dropdown', 'options'),
#      dash.Output('qc_scatter1_dropdown', 'options'),
#      dash.Output('qc_scatter2_dropdown', 'options'),
#      dash.Output('qc_scatter3_dropdown', 'options'),
#      dash.Output('qc_per_condition_dropdown2', 'options'),
#     [
#      dash.Input('add-row-button', 'n_clicks'),
#      dash.Input('qc_pattern_table', 'data'),
#      dash.Input('qc_pattern_table', 'active_cell'),
#     ],
#     [
#      dash.State('qc_pattern_table', 'columns')
#     ]
# )
# def update_qc_pattern_table(n_clicks, table_patterns, _, columns):

#     # print(n_clicks," ",config.qc_n_clicks_old)
#     if n_clicks > config.qc_n_clicks_old:
#         config.qc_n_clicks_old += 1
#         table_patterns.append({col['id']: '' for col in columns})

#     f_qc_update_patterns(config.adata, table_patterns)
#     table_patterns = f_qc_table_pattern(config.adata) #get the updated table
#     table_threshold = f_qc_table_threshold(config.adata) #update thresholds table
#     options = f_options(config.adata,"qc_")

#     return table_patterns, table_threshold, options, options, options, options, options

# @app.callback(
#      dash.Output('qc_thresholds_table', 'data'),
#     [
#         dash.Input('qc_pattern_table', 'data'),
#         dash.Input('qc_thresholds_table', 'data')        
#     ]
# )
# def update_qc_threshold_table(current_rows, table):

#     data = []
#     for m in table:
#         if m["Control measure"] in config.adata.uns["qc"].keys():
#             try:
#                 config.adata.uns["qc"][m["Control measure"]]["Minimum threshold"] = float(m["Minimum Threshold"])
#             except:
#                 config.adata.uns["qc"][m["Control measure"]]["Minimum threshold"] = config.adata.obs[m["Control measure"]].min()

#             try:
#                 config.adata.uns["qc"][m["Control measure"]]["Maximum threshold"] = float(m["Maximum Threshold"])
#             except:
#                 config.adata.uns["qc"][m["Control measure"]]["Maximum threshold"] = config.adata.obs[m["Control measure"]].max()

#     for i in config.adata.uns["qc"]:
#         data.append({'Control measure':i, 'Minimum Threshold':config.adata.uns["qc"][i]["Minimum threshold"], 'Maximum Threshold':config.adata.uns["qc"][i]["Maximum threshold"]})

#     return data

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
     dash.Output('data', 'children'),
    [
        dash.Input('doublets-button', 'n_clicks'),
    ]
)
def update_qc_threshold_per_condition_table(condition):
    if condition % 2 == 0:
        return [dbc.Button(id='doublets-button', n_clicks=0, children="Add Doublet Analysis",
                            size="lg",
                            style={
                                "background-color":"gray",
                                # 'width': '680px', 
                            }      
                            )]
    else:

        l = make_arguments("doublet_",doublet_args,title="Scrublet arguments")

        return  [
                    dbc.Row(
                        [
                            dbc.Col(
                                html.H1("Doublet Analysis")
                            ),
                            dbc.Col(
                                dbc.Button(id='doublets-button', n_clicks=1, children="Remove Doublet Analysis",
                                        size="lg",
                                        style={
                                            "background-color":"#343A40",
                                            'width': '280px', 
                                        }      
                                        ),
                                width=2
                            ),
                        ],
                        style={"margin-bottom":"1cm"}
                    ),
                    dbc.Row(
                            dcc.Dropdown(id='qc_doublet_scrolldown'),
                            style={"margin-bottom":"1cm"}
                        ),  
                ] + \
                [
                    dbc.Col(l,width=4,
                            style={"background-color":"lightgray"}
                            ),
                    dbc.Col(
                        dcc.Graph(id='qc_plot_doublet'),
                    ),  
                ]