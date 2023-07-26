import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import os
import scanpy as sc
from . import config
import pandas as pd
from .functions import *
import dash_table
import dash_renderjson

from app import app

# needed only if running this as a single page app
#external_stylesheets = [dbc.themes.LUX]

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# change to app.layout if running as single page app instead
def layout():

    config.h5ad_files = np.sort([f for f in os.listdir(config.folder_path) if f.endswith('.h5ad')])

    return dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(html.H1("Load Dataset"), width="auto"),
            ],
            justify="center",
            className="mb-4"
        ),
        dbc.Row(
            [
                dbc.Col(
                    html.Div(
                        dcc.Dropdown(
                            id='h5ad-dropdown',
                            options=[{'label': file, 'value': file} for file in config.h5ad_files],
                            value=config.old_selected_file,
                            clearable=False
                        )
                    )
                )
            ],
            justify="center",
            style={"margin-bottom":"1cm"}
        ),
        dbc.Row(
            [
                dbc.Col(html.H1(".X"), width="auto"),
            ],
            justify="left",
            className="mb-4"
        ),
        dbc.Row(
            [
                dash_table.DataTable(
                    id='table_x',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["cells","genes","dtype"]
                    ],
                    data=[{"cells": str(config.adata.X.shape[0]), "genes": str(config.adata.X.shape[1]), "dtype": str(type(config.adata.X))}],
                    editable=False,
                    row_deletable=False,
                    fixed_rows={'headers': True},
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                ),
            ],
            style={"margin-bottom":"1cm"}
        ),
        dbc.Row(
            [
                dbc.Col(html.H1(".obs"), width="auto"),
            ],
            justify="left",
            className="mb-4"
        ),
        dbc.Row(
            [
                dash_table.DataTable(
                    id='table_obs',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in config.adata.obs.columns.values
                    ],
                    data=config.adata.obs.to_dict("records"),
                    editable=True,
                    row_deletable=False,
                    fixed_rows={'headers': True},
                    style_table={'overflowY': 'auto', 'height': '500px', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.H1(".var"), width="auto"),
            ],
            justify="left",
            className="mb-4"
        ),
        dbc.Row(
            [
                dash_table.DataTable(
                    id='table_var',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in config.adata.var.columns.values
                    ],
                    data=config.adata.var.sort_values("Gene").to_dict("records"),
                    editable=True,
                    row_deletable=False,
                    fixed_rows={'headers': True},
                    style_table={'overflowY': 'auto', 'height': '500px', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                ),
            ]
        ),
        dbc.Row(
            [
                dbc.Col(html.P("Gene lists."), width="auto"),
            ],
            justify="left",
            className="mb-4"
        ),
        dbc.Row(
        dash_table.DataTable(
                id='qc_pattern_table',
                columns=[
                    {"name": i, "id": i, "deletable": False, "editable": True} for i in config.qc_df_patterns.columns
                ],
                data=f_qc_table_pattern(config.adata),
                editable=True,
                row_deletable=True,
                # row_selectable="multi",
                style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
            ),
        ),
        dbc.Row(
            dbc.Button('Add Row', id='add-row-button', n_clicks=0,
                       style={
                           "background-color":"#343A40",
                        }
                        ),
            style={"margin-bottom":"1cm"}
        ),
        dbc.Row(
            [
                dbc.Col(html.H1(".obsm"), width="auto"),
            ],
            justify="left",
            className="mb-4"
        ),
        dbc.Row(
            [
                dash_table.DataTable(
                    id='table_obsm',
                    columns=[
                        {"name": i, "id": i, "deletable": False, "editable": False} for i in ["name","cells","vars","dtype"]
                    ],
                    data=[{"name": i, "cells": str(config.adata.obsm[i].shape[0]), "vars": str(config.adata.obsm[i].shape[1]), "dtype": str(type(config.adata.obsm[i]))}
                          for i in config.adata.obsm
                          ],
                    editable=False,
                    row_deletable=False,
                    fixed_rows={'headers': True},
                    style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                    style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto', 'minWidth': 90}
                ),
            ],
            style={"margin-bottom":"1cm"}
        ),
        dbc.Row(
            [
                dbc.Col(html.H1(".uns"), width="auto"),
            ],
            justify="left",
            className="mb-4"
        ),
        dbc.Row(
            dash_renderjson.DashRenderjson(id="output_uns", data=json_serializable(config.adata.uns), max_depth=1, invert_theme=True)
        ),
        dbc.Row(
            [
                dbc.Button(id='gene-list-save-button', n_clicks=0, children="Save",
                            size="lg",
                            style={
                                "background-color":"#343A40",
                                   'width': '280px', 
                            }      
                            )
            ],
            justify="center",
            className="mb-4"
        ),

    ]
)

# @app.callback(
#     dash.dependencies.Output('h5ad-dropdown', 'clearable'),
#     [dash.dependencies.Input('h5ad-dropdown', 'value')],
# )
# def load_data(selected_file):
    
#     # Upload the file if changed
#     if selected_file != config.old_selected_file:
#         config.file_path = os.path.join(config.folder_path, selected_file)
#         config.old_selected_file = selected_file
#         config.adata = sc.read(config.file_path)

#         config.f_qc(config.adata)

#     return False

@app.callback(
     dash.Output('add-row-button', 'n_clicks'),
     dash.Output('qc_pattern_table', 'data'),
     dash.Output('table_x', 'data'),
     dash.Output('table_obs', 'columns'),
     dash.Output('table_obs', 'data'),
     dash.Output('table_var', 'columns'),
     dash.Output('table_var', 'data'),
     dash.Output('output_uns', 'data'),
    [
     dash.Input('h5ad-dropdown', 'value'),
     dash.Input('add-row-button', 'n_clicks'),
     dash.Input('qc_pattern_table', 'data'),
     dash.Input('qc_pattern_table', 'active_cell'),
    ],
    [
     dash.State('qc_pattern_table', 'columns')
    ]
)
def update_pattern_table(selected_file, n_clicks, table_patterns, _, columns):

    if selected_file != config.old_selected_file:
        config.file_path = os.path.join(config.folder_path, selected_file)
        config.old_selected_file = selected_file
        config.adata = sc.read(config.file_path)

        config.qc_n_clicks_old = np.Inf

        config.f_qc(config.adata)

    add_clicks = 0

    if n_clicks == 0:
        config.qc_n_clicks_old = 1
        add_clicks = 1
    elif n_clicks < config.qc_n_clicks_old:
        config.qc_n_clicks_old = 1
        add_clicks = 1
    elif n_clicks > config.qc_n_clicks_old:
        config.qc_n_clicks_old += 1
        add_clicks = n_clicks
        table_patterns.append({col['id']: '' for col in columns})
        f_update_patterns(config.adata, table_patterns)
    else:
        add_clicks = n_clicks
        f_update_patterns(config.adata, table_patterns)
    
    table_patterns = f_qc_table_pattern(config.adata) #get the updated table

    columns_obs=[
                {"name": i, "id": i, "deletable": False, "editable": False} for i in config.adata.obs.columns.values
            ]
    data_obs=config.adata.obs.to_dict("records")

    columns_var=[
                {"name": i, "id": i, "deletable": False, "editable": False} for i in config.adata.var.columns.values
            ]
    data_var=config.adata.var.sort_values("Gene").to_dict("records")

    x = [{"cells": str(config.adata.X.shape[0]), "genes": str(config.adata.X.shape[1]), "dtype": str(type(config.adata.X))}]

    return add_clicks, table_patterns, x, columns_obs, data_obs, columns_var, data_var, json_serializable(config.adata.uns)

@app.callback(
     dash.Output('gene-list-save-button', 'children'),
    [
     dash.Input('gene-list-save-button', 'n_clicks'),
    ]
)
def save_gene_lists(n_clicks):

    config.adata.write(config.file_path)

    return "Save"