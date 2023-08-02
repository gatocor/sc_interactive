import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash as dcc
import os
import scanpy as sc
from . import config
import pandas as pd
from .functions import *

from app import app

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

error_modal = dbc.Modal(
    [
        dbc.ModalHeader("Error"),
        dbc.ModalBody(id="error-message",
                      children=dcc.Dropdown(id="gene-dropdown",
                                            value=None,
                                            options=[]#config.adata.var.columns.values
                                            )
                    ),
        dbc.ModalFooter(dbc.Button("Save", id="error-close", className="ml-auto")),
    ],
    id="error-modal",
    size="sm",
)

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
                            dcc.Dropdown(
                                id='h5ad-dropdown',
                                options=[{'label': file, 'value': file} for file in config.h5ad_files],
                                value=config.selected_file,
                                clearable=False
                            )
                    ),
                    dbc.Col(
                            dbc.Button(
                                id='h5ad-load-button', n_clicks=0, children="Load",
                                size="lg",
                                style={
                                    "background-color":"#343A40",
                                    'width': '280px', 
                                }                        )
                    )
                ],
                justify="center",
                style={"margin-bottom":"1cm"}
            ),
            error_modal,
            dbc.Row(id="adata_info",children=make_adata_layout(config.adata))
        ],
        fluid=True
    )

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
     dash.Input('add-row-button', 'n_clicks'),
     dash.Input('qc_pattern_table', 'data'),
     dash.Input('qc_pattern_table', 'active_cell'),
    ],
    [
     dash.State('h5ad-dropdown', 'value'),
     dash.State('qc_pattern_table', 'columns')
    ],
)
def update_pattern_table(n_clicks, table_patterns, _, selected_file, columns):

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
    data_var=config.adata.var.sort_values(config.adata.uns["GeneNamesKey"]).to_dict("records")

    x = [{"cells": str(config.adata.X.shape[0]), "genes": str(config.adata.X.shape[1]), "dtype": str(type(config.adata.X))}]

    return add_clicks, table_patterns, x, columns_obs, data_obs, columns_var, data_var, json_serializable(config.adata.uns)

@app.callback(
     dash.Output('adata_info', 'children'),
     dash.Output('error-modal', 'is_open'),
     dash.Output('error-message', 'children'),
     dash.Output('error-close', 'n_clicks'),
    [
     dash.Input('h5ad-load-button', 'n_clicks'),
     dash.Input('error-close', 'n_clicks'),
    ],
    [
     dash.State('gene-dropdown', 'value'),
     dash.State('h5ad-dropdown', 'value'),
    ],
    prevent_initial_call=True
)
def load_adata_info(n, n_clicks, gene, selected_file):

    if selected_file != None:

        config.file_path = os.path.join(config.folder_path, selected_file)
        config.selected_file = selected_file
        config.adata = sc.read(config.file_path)
        config.f_qc_base(config.adata)

        if "GeneNamesKey" in config.adata.uns.keys():

            None

        elif type(gene) != type(None):

            config.adata.uns["GeneNamesKey"] = gene
            gene = None

        elif ("GeneNamesKey" not in config.adata.uns) or (config.adata.uns["GeneNamesKey"] == None):

            layout = [
                html.Div("Your dataset does not have specified which key in .var is the gene names key. Please specify one before continuing."),
                        dcc.Dropdown(id="gene-dropdown",
                                    value=None,
                                    options=config.adata.var.columns.values
                                    )
            ]
            return [], True, layout, n_clicks

        layout = make_adata_layout(config.adata)

        return layout, False, "", n_clicks
    else:
        return [], False, "", 0


@app.callback(
     dash.Output('gene-list-save-button', 'children'),
    [
     dash.Input('gene-list-save-button', 'n_clicks'),
    ],
    prevent_initial_call=True
)
def save_gene_lists(n_clicks):

    config.adata.write(config.file_path)

    return "Save"