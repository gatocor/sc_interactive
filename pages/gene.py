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

layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(html.H1("Additional tags"), width="auto"),
            ],
            justify="center",
            className="mb-4"
        ),
        dbc.Row(
            [
                dbc.Col(html.P("In here you can put any lists."), width="auto"),
            ],
            justify="left",
            className="mb-4"
        ),
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
     dash.Output('qc_pattern_table', 'data'),
    [
     dash.Input('add-row-button', 'n_clicks'),
     dash.Input('qc_pattern_table', 'data'),
     dash.Input('qc_pattern_table', 'active_cell'),
    ],
    [
     dash.State('qc_pattern_table', 'columns')
    ]
)
def update_pattern_table(n_clicks, table_patterns, _, columns):

    # print(n_clicks," ",config.qc_n_clicks_old)
    if n_clicks > config.qc_n_clicks_old:
        config.qc_n_clicks_old += 1
        table_patterns.append({col['id']: '' for col in columns})

    print("Hola ", n_clicks, config.adata.uns["gene_lists"], table_patterns)

    if n_clicks != 0:
        f_update_patterns(config.adata, table_patterns)
    elif n_clicks < config.qc_n_clicks_old:
        config.qc_n_clicks_old = 0
    
    table_patterns = f_qc_table_pattern(config.adata) #get the updated table

    print("Hola2 ",table_patterns)

    return table_patterns
