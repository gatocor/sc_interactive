import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import os
import scanpy as sc
from . import config
import pandas as pd
from .functions import *

from app import app

# needed only if running this as a single page app
#external_stylesheets = [dbc.themes.LUX]

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)

# change to app.layout if running as single page app instead
layout = dbc.Container(
    [
        dbc.Row(
            [
                dbc.Col(html.H1("Home"), width="auto"),
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
        ),
    ]
)

@app.callback(
    dash.dependencies.Output('h5ad-dropdown', 'clearable'),
    [dash.dependencies.Input('h5ad-dropdown', 'value')],
)
def load_data(selected_file):
    
    # Upload the file if changed
    if selected_file != config.old_selected_file:
        file_path = os.path.join(config.folder_path, selected_file)
        config.old_selected_file = selected_file
        config.adata = sc.read(file_path)

        config.f_qc(config.adata)

    return False