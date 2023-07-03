import os
import dash
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import scanpy as sc
import plotly.graph_objs as go
import numpy as np
import pandas as pd
from . import config

from app import app

#Lists
control_options = ["n_genes_by_counts","log1p_n_genes_by_counts","total_counts","log1p_total_counts"]

data = {
    'Name': ['Alice', 'Bob', 'Charlie'],
    'Age': [25, 30, 35],
    'City': ['New York', 'London', 'Paris']
}
df = pd.DataFrame(data)

layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(html.H1("Quality control"), width="auto"),
            ],
            justify="center",
            className="mb-4"
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
                                        id='control_plot_type_dropdown',
                                        options=["histogram","scatter","violin"],
                                        value="histogram",
                                        placeholder="Select a column",
                                        clearable=False
                                    )
                                ),
                            ]
                        ),
                    ]),
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
                            html.Div(
                                dcc.Dropdown(
                                    id='control_var1_dropdown',
                                    options=control_options,
                                    value="total_counts",
                                    placeholder="Select a column",
                                    clearable=False
                                )
                            ),
                        ),
                        dbc.Col(
                            html.Div(
                                dcc.Input(
                                    id='var1_min-threshold-input',
                                    type='number',
                                    placeholder='mt',
                                    value=0,
                                    debounce=True
                                )
                            ),
                        ),
                        dbc.Col(
                            html.Div(
                                dcc.Input(
                                    id='var1_max-threshold-input',
                                    type='number',
                                    placeholder='mt',
                                    value=0,
                                    debounce=True
                                )
                            ),
                        )
                    ]),
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
                            html.Div(
                                dcc.Dropdown(
                                    id='control_var2_dropdown',
                                    options=control_options,
                                    value="total_counts",
                                    placeholder="Select a column",
                                    clearable=False,
                                    disabled=True,
                                )
                            )
                        ),
                        dbc.Col(
                            html.Div(
                                dcc.Input(
                                    id='var2_min-threshold-input',
                                    type='number',
                                    placeholder='mt',
                                    value=0,
                                    debounce=True,
                                    disabled=True,
                                )
                            ),
                        ),
                        dbc.Col(
                            html.Div(
                                dcc.Input(
                                    id='var2_max-threshold-input',
                                    type='number',
                                    placeholder='mt',
                                    value=0,
                                    debounce=True,
                                    disabled=True,
                                )
                            ),
                        )
                    ]),
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
                            html.Div(
                                dcc.Dropdown(
                                    id='control_var3_dropdown',
                                    options=control_options,
                                    value="total_counts",
                                    placeholder="Select a column",
                                    clearable=False,
                                    disabled=True,
                                )
                            ),
                        )
                    ]),
                    dcc.Graph(
                        id='table',
                        figure=go.Figure(
                            data=[go.Table(
                                header=dict(values=df.columns),
                                cells=dict(values=[df[col] for col in df.columns])
                            )],
                            layout=go.Layout(
                                title='Data Table'
                            )
                        )
                    )
                ],
                width=6
                ),
                dbc.Col(
                    dcc.Graph(id='histogram-plot'),
                    width=6
                ),
            ],
            justify="center",
            align="top"
        ),
    ],
    fluid=True
)

# Step 3: Add a callback function to update the column dropdown options based on the selected .h5ad file
# @app.callback(
#     dash.dependencies.Output('column-dropdown', 'options'),
#     [dash.dependencies.Input('h5ad-dropdown', 'value'),
#      dash.dependencies.Input('column-dropdown', 'value')]
# )
# def update_column_dropdown(selected_file, var1_select):
#     global adata  # Access the global variable
#     column_names = adata.obs.columns.tolist()

#     return [{'label': column, 'value': column} for column in column_names]

# Step 4: Add a callback function to update the histogram plot based on the selected column
@app.callback(
    dash.dependencies.Output('histogram-plot', 'figure'),
    dash.dependencies.Output('control_var2_dropdown', 'disabled'),
    dash.dependencies.Output('var2_min-threshold-input', 'disabled'),
    dash.dependencies.Output('var2_max-threshold-input', 'disabled'),
    dash.dependencies.Output('control_var3_dropdown', 'disabled'),
    [#dash.dependencies.Input('h5ad-dropdown', 'value'),
     dash.dependencies.Input('control_plot_type_dropdown', 'value'),
     dash.dependencies.Input('control_var1_dropdown', 'value'),
     dash.dependencies.Input('var1_min-threshold-input', 'value'),
     dash.dependencies.Input('var1_max-threshold-input', 'value'),  
     dash.dependencies.Input('control_var2_dropdown', 'value'),
     dash.dependencies.Input('var2_min-threshold-input', 'value'),
     dash.dependencies.Input('var2_max-threshold-input', 'value'),     
     dash.dependencies.Input('control_var3_dropdown', 'value'),
     ],
)
def update_qc(plot_type, var1_select, var1_min_threshold, var1_max_threshold, var2_select, var2_min_threshold, var2_max_threshold, var3_select):
    
    # Perform necessary computations or filtering on the adata object
    if "total_counts" not in config.adata.obs.columns.values:
        sc.pp.calculate_qc_metrics(config.adata,percent_top=(3,),inplace=True)#np.round(np.int,np.linspace(1,adata.shape[1],5)))

    # Example: Plot histogram of the selected column from adata.obs
    var1_selected_data = config.adata.obs[var1_select].values
    var2_selected_data = config.adata.obs[var2_select].values
    var3_selected_data = config.adata.obs[var3_select].values

    #Plot_type
    control_var2_dropdown = True
    var2_var1_min_threshold_input_disabled = True
    var2_var1_max_threshold_input_disabled = True
    control_var3_dropdown = True
    if plot_type == "histogram":
        control_var2_dropdown = True
        var2_var1_min_threshold_input_disabled = True
        var2_var1_max_threshold_input_disabled = True
        control_var3_dropdown = True

        # Create a vertical line at the specified input value
        hist_values, hist_bins = np.histogram(var1_selected_data, bins=30)
        tallest_bin_height = np.max(hist_values)
        var1_vertical_line_min = go.Scatter(
            x=[var1_min_threshold, var1_min_threshold],
            y=[0, tallest_bin_height],
            mode='lines',
            name='Min threshold',
            line=dict(color='red', dash='dash')
        )
        var1_vertical_line_max = go.Scatter(
            x=[var1_max_threshold, var1_max_threshold],
            y=[0, tallest_bin_height],
            mode='lines',
            name='Max threshold',
            line=dict(color='green', dash='dash')
        )

        data = [
            go.Histogram(
                x=var1_selected_data,
                nbinsx=30,
                name='Histogram',
                marker=dict(color='blue'),
                opacity=0.7
            ),
            var1_vertical_line_min,
            var1_vertical_line_max
        ]

        layout = {
            'title': f'Histogram of {var1_select}',
            'xaxis': {'title': var1_select},
            'yaxis': {'title': 'Count'},
            'barmode': 'overlay',
            'width':900,
            'height':800,
        }

    elif plot_type == "scatter":
        control_var2_dropdown = False
        var2_var1_min_threshold_input_disabled = False
        var2_var1_max_threshold_input_disabled = False
        control_var3_dropdown = False

        # Create a vertical line at the specified input value
        var1_vertical_line_min = go.Scatter(
            x=[var1_min_threshold, var1_min_threshold],
            y=[np.min(var2_selected_data), np.max(var2_selected_data)],
            mode='lines',
            name='Min threshold',
            line=dict(color='red', dash='dash')
        )
        var1_vertical_line_max = go.Scatter(
            x=[var1_max_threshold, var1_max_threshold],
            y=[np.min(var2_selected_data), np.max(var2_selected_data)],
            mode='lines',
            name='Max threshold',
            line=dict(color='green', dash='dash')
        )

        # Create a vertical line at the specified input value
        var2_vertical_line_min = go.Scatter(
            y=[var2_min_threshold, var2_min_threshold],
            x=[np.min(var1_selected_data), np.max(var1_selected_data)],
            mode='lines',
            name='Min threshold',
            line=dict(color='red', dash='dot')
        )
        var2_vertical_line_max = go.Scatter(
            y=[var2_max_threshold, var2_max_threshold],
            x=[np.min(var1_selected_data), np.max(var1_selected_data)],
            mode='lines',
            name='Max threshold',
            line=dict(color='green', dash='dot')
        )

        data = [
            go.Scatter(
                x=var1_selected_data,
                y=var2_selected_data,
                marker=dict(
                    color=var3_selected_data,
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
            'title': f'Scatter of {var1_select} vs {var2_select}',
            'xaxis': {'title': var1_select},
            'yaxis': {'title': var2_select,
                        # 'scaleanchor': 'x',
                        'scaleratio': 1
                      },
            'barmode': 'overlay',
            'width':900,
            'height':800,
        }
        
    elif plot_type == "violin":
        control_var2_dropdown = False
        var2_var1_min_threshold_input_disabled = False
        var2_var1_max_threshold_input_disabled = False
        control_var3_dropdown = True

        data = [
            go.Violin(x=var3_selected_data, y=var1_select), # Include hover data
            var1_vertical_line_min,
            var1_vertical_line_max
        ]

    return {
        'data': data,
        'layout': layout,
    }, control_var2_dropdown, var2_var1_min_threshold_input_disabled, var2_var1_max_threshold_input_disabled, control_var3_dropdown

# # Step 5: Run the app
# if __name__ == '__main__':
#     app.run_server(debug=True)
