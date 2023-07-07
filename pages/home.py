import dash
import dash_html_components as html
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import os
import dash_table
import scanpy as sc
from . import config
import pandas as pd

from app import app

# needed only if running this as a single page app
#external_stylesheets = [dbc.themes.LUX]

#app = dash.Dash(__name__, external_stylesheets=external_stylesheets)
data = {
    'Pattern': ['mt-'],
    'Genes following the pattern': [str([0])]
}
df = pd.DataFrame(data)

n_clicks_old = 0

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
        dash_table.DataTable(
                id='table',
                columns=[
                    {"name": i, "id": i, "deletable": False, "editable": True if i == 'Pattern' else False} for i in df.columns
                ],
                data=df.to_dict('records'),
                editable=True,
                row_deletable=True,
                style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
            ),
        html.Button('Add Row', id='add-row-button', n_clicks=0),
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

    return False

# Callback to handle adding and deleting rows
@app.callback(
     dash.Output('table', 'data'),
    [dash.Input('add-row-button', 'n_clicks'),
     dash.Input('table', 'data')],
    [dash.State('table', 'data_previous'),
     dash.State('table', 'data'),
     dash.State('table', 'columns')]
)
def update_table(n_clicks, data_, previous_rows, current_rows, columns):
    global n_clicks_old
    if n_clicks > n_clicks_old:
        current_rows.append({col['id']: '' for col in columns})
        n_clicks_old += 1

    if previous_rows is not None:
        deleted_rows = [row for row in previous_rows if row not in current_rows]
        if deleted_rows:
            print(f"Deleted rows: {deleted_rows}")

    # print(current_rows)
    for i in range(len(current_rows)):
        pattern = current_rows[i]['Pattern']
        if pattern != '':
            c = [j for j in config.adata.var.loc[:,"Gene"] if pattern in j]
        else:
            c = ''
        current_rows[i] = {'Pattern':pattern, 'Genes following the pattern':str(c)}

    return current_rows
