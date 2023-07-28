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
            children = [],
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
                                    options=[i for i in config.adata.uns["qc"]],
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
    dash.Output('table_qc_metrics','data'),
    dash.Output('qc_thresholds_table', 'data'),
    dash.Output("global_qc_plots","children"),
    dash.Output("dropdown_add_metric","options"),
    dash.Output("dropdown_add_metric","value"),
    dash.Output('qc_scatter1_dropdown', 'options'),
    dash.Output('qc_scatter2_dropdown', 'options'),
    dash.Output('qc_scatter3_dropdown', 'options'),
    [
        dash.Input('add_metric', 'n_clicks'),
        dash.Input('table_qc_metrics', 'data'),
        dash.Input('qc_thresholds_table', 'data'),
    ],
    [
        dash.State('dropdown_add_metric', 'value'),
    ]
)
def update_qc_global_histogram(n_clicks, qc_metrics, data, value):
    
    l = []

    if value != None:
        qc_metrics.append({"Name":str(value)})

    f_qc(config.adata, qc_metrics, data)

    for var_selected_data in [i for i in config.adata.uns["qc"]]:
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

        l += [
                dbc.Row(
                    [dcc.Graph(id="Holi",
                            figure={
                                    "data":[
                                        go.Histogram(
                                            x=config.adata.obs[var_selected_data].values,
                                            nbinsx=30,
                                            name='Histogram',
                                            marker=dict(color='blue'),
                                            opacity=0.7
                                        ),
                                        var1_vertical_line_min,
                                        var1_vertical_line_max
                                    ],
                                    "layout":{
                                            'title': f'Histogram of {var_selected_data}',
                                            'xaxis': {'title': var_selected_data},
                                            'yaxis': {'title': 'Count'},
                                            'barmode': 'overlay',
                                            'width':1500,
                                            'height':400,
                                    }
                                }
                    )],
                    justify="center",
                    style={'width': '90%', 'margin': 'auto'}
                )
            ]
        
    options = [str(i) for i in config.adata.var.columns.values if (config.adata.var.dtypes[i] in [bool]) and (i not in config.adata.uns["qc"])]
        
    options_dropdown = list(config.adata.uns["qc"].keys())

    return f_qc_table_metrics(config.adata), f_qc_table_threshold(config.adata), l, options, None, options_dropdown, options_dropdown, options_dropdown #{'data': data,'layout': layout}

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
    ]
)
def make_qc_per_condition_(condition):
        
    return make_qc_per_condition(config.adata)

@app.callback(
     dash.Output('per_condition_plot', 'children'),
     dash.Output('table_qc_per_condition_metrics', 'data'),
     dash.Output('dropdown_add_per_condition_metrics', 'value'),
    [
        dash.Input('add_qc_per_condition-button', 'n_clicks'),
        dash.Input('table_qc_per_condition_metrics', 'data'),
    ],
    dash.State('dropdown_add_per_condition_metrics', 'value'),
)
def update_qc_threshold_per_condition_table(condition, data, value):

    conditions = [i["Name"] for i in data]
    if value != None and value not in conditions:
        data.append({"Name":str(value)})

    conditions = [i["Name"] for i in data]

    l = []

    for var_selected_data in [i for i in config.adata.uns["qc"]]:          
        
        if len(data) != 0:                
            for condition in conditions:
                if condition not in config.adata.uns["qc"][var_selected_data]["per_condition"].keys():
                    config.adata.uns["qc"][var_selected_data]["per_condition"][condition] = \
                        {
                            str(ci):{
                                "Min":config.adata.uns["qc"][var_selected_data]["Minimum threshold"],
                                "Max":config.adata.uns["qc"][var_selected_data]["Maximum threshold"]                                                              
                            }
                            for ci in np.unique(config.adata.obs[condition].values)
                        }

    for condition in conditions:
        l += [
                dbc.Row(
                    html.H1(condition),
                    justify="left"
                )
        ]
        for var_selected_data in [i for i in config.adata.uns["qc"]]:                
            #Plot_type
            # Create a vertical line at the specified input value
            var1_vertical_line_min = go.Scatter(
                x=config.adata.obs[condition].values,
                y=[config.adata.uns["qc"][var_selected_data]["Minimum threshold"] for i in range(len(config.adata.obs[var_selected_data].values))],
                mode='lines',
                name='Min threshold',
                line=dict(color='red')
            )
            var1_vertical_line_max = go.Scatter(
                x=config.adata.obs[condition].values,
                y=[config.adata.uns["qc"][var_selected_data]["Maximum threshold"] for i in range(len(config.adata.obs[var_selected_data].values))],
                mode='lines',
                name='Max threshold',
                line=dict(color='green')
            )

            l += [
                    dbc.Col(
                        [
                            dcc.Graph(id="Holi",
                                figure={
                                        "data":[
                                            go.Violin(
                                                x=config.adata.obs[condition].values,
                                                y=config.adata.obs[var_selected_data].values,
                                                name='Violin',
                                                # marker=dict(color='blue'),
                                                opacity=0.7
                                            ),
                                            var1_vertical_line_min,
                                            var1_vertical_line_max
                                        ],
                                        "layout":{
                                                'title': f'{var_selected_data}',
                                                'xaxis': {'title': condition},
                                                'yaxis': {'title': 'Count'},
                                                'barmode': 'overlay',
                                                'width':1500,
                                                'height':400,
                                        }
                                    }
                            )
                        ],
                        # justify="center",
                        # style={'width': '90%', 'margin': 'auto'}
                    )
                ]
            
            l += [
                dbc.Col(
                    dash_table.DataTable(
                            id='table_qc_per_condition_'+condition,
                            columns=[
                                {"name": str(i), "id": str(i), "deletable": False, "editable": j} for i,j in zip(["Condition","Min","Max"],[False,True,True])
                            ],
                            data=[{"Condition":i,"Min":"i","Max":i} for i in np.unique(config.adata.obs[condition].values)],
                            editable=True,
                            row_deletable=False,
                            style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                            style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                        ),
                )
            ]

    return l, data, None

@app.callback(
     dash.Output('doublets', 'children'),
    [
        dash.Input('doublets-button', 'n_clicks'),
    ]
)
def make_qc_doublets(condition):
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