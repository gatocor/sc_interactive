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

from app import app

#Lists
layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(html.H1("Additional tags"), width="auto"),
            ],
            justify="center",
            className="mb-4"
        ),
        dash_table.DataTable(
                id='qc_pattern_table',
                columns=[
                    {"name": i, "id": i, "deletable": False, "editable": True if i == 'Pattern' else False} for i in config.qc_df_patterns.columns
                ],
                data=config.qc_df_patterns.to_dict('records'),
                editable=True,
                row_deletable=True,
                style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
            ),
        dbc.Row(
            dbc.Button('Add Row', id='add-row-button', n_clicks=0,
                       style={
                           "background-color":"#343A40",
                        }
                        ),
        ),
        dbc.Row(
            [
                dbc.Col(html.H1("Global Quality control"), width="auto"),
            ],
            justify="center",
            className="mb-4"
        ),
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
                                        options=[],
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
                                            options=[],
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
                                            options=[],
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
                                            options=[],
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
                value=[i for i in config.adata.obs.columns.values if i.startswith("condition_")][0],
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
@app.callback(
     dash.Output('qc_pattern_table', 'data'),
     dash.Output('qc_hist_dropdown', 'options'),
     dash.Output('qc_scatter1_dropdown', 'options'),
     dash.Output('qc_scatter2_dropdown', 'options'),
     dash.Output('qc_scatter3_dropdown', 'options'),
     dash.Output('qc_per_condition_dropdown2', 'options'),
    [dash.Input('add-row-button', 'n_clicks'),
     dash.Input('qc_pattern_table', 'data')],
    [dash.State('qc_pattern_table', 'data_previous'),
     dash.State('qc_pattern_table', 'data'),
     dash.State('qc_pattern_table', 'columns')]
)
def update_qc_pattern_table(n_clicks, data_, previous_rows, current_rows, columns):
    # print(current_rows)
    if n_clicks > config.qc_n_clicks_old:
        current_rows.append({col['id']: '' for col in columns})
        config.qc_n_clicks_old += 1

    if previous_rows is not None:
        deleted_rows = [row for row in previous_rows if row not in current_rows]
        if deleted_rows:
            print(f"Deleted rows: {deleted_rows}")

    # print(current_rows)
    for i in range(len(current_rows)):
        pattern = current_rows[i]['Pattern']
        if pattern != '':
            c = [j for j in config.adata.var.loc[:,"Gene"] if j.startswith(pattern)]
        else:
            c = ''
        current_rows[i] = {'Pattern':pattern, 'Genes following the pattern':str(c)}
        
        f_qc_pattern(config.adata,pattern)

    # print(config.adata.obs.columns.values)
    options = [i for i in config.adata.obs.columns.values if i.startswith("qc_")]

    return current_rows, options, options, options, options, options

@app.callback(
     dash.Output('qc_thresholds_table', 'data'),
    [
        dash.Input('qc_pattern_table', 'data'),
        dash.Input('qc_thresholds_table', 'data')        
    ]
)
def update_qc_threshold_table(current_rows, table):

    data = []
    for m in table:
        if m["Control measure"] in config.adata.uns["qc"].keys():
            try:
                config.adata.uns["qc"][m["Control measure"]]["Minimum threshold"] = float(m["Minimum Threshold"])
            except:
                config.adata.uns["qc"][m["Control measure"]]["Minimum threshold"] = config.adata.obs[m["Control measure"]].min()

            try:
                config.adata.uns["qc"][m["Control measure"]]["Maximum threshold"] = float(m["Maximum Threshold"])
            except:
                config.adata.uns["qc"][m["Control measure"]]["Maximum threshold"] = config.adata.obs[m["Control measure"]].max()

    for i in config.adata.uns["qc"]:
        data.append({'Control measure':i, 'Minimum Threshold':config.adata.uns["qc"][i]["Minimum threshold"], 'Maximum Threshold':config.adata.uns["qc"][i]["Maximum threshold"]})

    return data

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
        return  [
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
                    ),
                    dbc.Row(
                        dcc.Dropdown(id='qc_doublet_scrolldown'),
                    ),  

                    dbc.Col(
                        [
                            dbc.Tooltip(
                                "float, optional (default: 1.0) Rate for sampling UMIs when creating synthetic doublets. If 1.0, each doublet is created by simply adding the UMIs from two randomly sampled observed transcriptomes. For values less than 1, the UMI counts are added and then randomly sampled at the specified rate.",
                                target="doublet_arg1",
                                placement="bottom",
                            ),
                            html.Label("synthetic_doublet_umi_subsampling",id="doublet_arg1"),
                            dbc.Input(id="doublet_arg1",value=1.0,type="number")
                        ]
                    ),

        # use_approx_neighbors : bool, optional (default: True)
        #     Use approximate nearest neighbor method (annoy) for the KNN 
        #     classifier.

        # distance_metric : str, optional (default: 'euclidean')
        #     Distance metric used when finding nearest neighbors. For list of
        #     valid values, see the documentation for annoy (if `use_approx_neighbors`
        #     is True) or sklearn.neighbors.NearestNeighbors (if `use_approx_neighbors`
        #     is False).
            
        # get_doublet_neighbor_parents : bool, optional (default: False)
        #     If True, return the parent transcriptomes that generated the 
        #     doublet neighbors of each observed transcriptome. This information can 
        #     be used to infer the cell states that generated a given 
        #     doublet state.

        # min_counts : float, optional (default: 3)
        #     Used for gene filtering prior to PCA. Genes expressed at fewer than 
        #     `min_counts` in fewer than `min_cells` (see below) are excluded.

        # min_cells : int, optional (default: 3)
        #     Used for gene filtering prior to PCA. Genes expressed at fewer than 
        #     `min_counts` (see above) in fewer than `min_cells` are excluded.

        # min_gene_variability_pctl : float, optional (default: 85.0)
        #     Used for gene filtering prior to PCA. Keep the most highly variable genes
        #     (in the top min_gene_variability_pctl percentile), as measured by 
        #     the v-statistic [Klein et al., Cell 2015].

        # log_transform : bool, optional (default: False)
        #     If True, log-transform the counts matrix (log10(1+TPM)). 
        #     `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
        #     reduction, unless `mean_center` is True.

        # mean_center : bool, optional (default: True)
        #     If True, center the data such that each gene has a mean of 0.
        #     `sklearn.decomposition.PCA` will be used for dimensionality
        #     reduction.

        # normalize_variance : bool, optional (default: True)
        #     If True, normalize the data such that each gene has a variance of 1.
        #     `sklearn.decomposition.TruncatedSVD` will be used for dimensionality
        #     reduction, unless `mean_center` is True.

        # n_prin_comps : int, optional (default: 30)
        #     Number of principal components used to embed the transcriptomes prior
        #     to k-nearest-neighbor graph construction.

        # svd_solver : str, optional (default: 'arpack')
        #     SVD solver to use. See available options for 
        #     `svd_solver` from `sklearn.decomposition.PCA` or
        #     `algorithm` from `sklearn.decomposition.TruncatedSVD`

                    dbc.Row(
                        dcc.Graph(id='qc_plot_doublet'),
                    ),  
                ]


