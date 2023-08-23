from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import os
import scanpy as sc

# must add this line in order for the app to be deployed successfully on Heroku
from app import server
from app import app
# import all pages in the app
from pages import config, home, graph

navbar = dbc.Navbar(
    dbc.Container(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/assets/scRNAseq.png", height="60px")),
                        dbc.Col(dbc.NavbarBrand("scRNAseq analysis", id="nav_brand", className="ml-2")),
                    ],
                    align="center",
                ),
                href="/home",
            ),
            dbc.Nav([
                dbc.DropdownMenu([
                    dbc.DropdownMenuItem("Dataset", href="/home"),                     
                    # dbc.DropdownMenuItem("Quality Control", href="/quality_control"),                  
                    dbc.DropdownMenuItem("Graphical analysis", href="/graph"),                     
                ],
                in_navbar = True,
                label = "Menu",
                color="dark"
                ),
            ])
        ]
    ),
    color="dark",
    dark=True,
    className="mb-4",
)

def toggle_navbar_collapse(n, is_open):
    if n:
        return not is_open
    return is_open

for i in [2]:
    app.callback(
        Output(f"navbar-collapse{i}", "is_open"),
        [Input(f"navbar-toggler{i}", "n_clicks")],
        [State(f"navbar-collapse{i}", "is_open")],
    )(toggle_navbar_collapse)

# embedding the navigation bar
app.layout = html.Div([
    dcc.Location(id='url', refresh=False),
    navbar,
    html.Div(id='page-content')
])

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/home':

        return home.layout()
    
    # elif pathname == '/quality_control':
        
    #     if config.CACHEFOLDER not in config.file_path:
    #         config.file_path = config.file_path.split(config.CACHEFOLDER)[-1]
    #         config.adata = sc.read(config.file_path)

    #     return quality_control.layout()
    elif pathname == '/graph':

        # if config.CACHEFOLDER not in config.file_path:
        #     config.file_path = config.CACHEFOLDER+config.file_path
        #     config.adata = sc.read(config.file_path)

        return graph.layout()

@app.callback(Output('nav_brand', 'children'),
              [Input('h5ad-load-button','n_clicks')],
              [State('h5ad-dropdown', 'value')])
def display_page(_,name):
    if type(name) != type(None):
        return "scRNAseq analysis: "+name.split(".h5ad")[0]
    else:
        return "scRNAseq analysis"

if __name__ == '__main__':
    app.run_server(debug=True)