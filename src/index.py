from dash import dcc
from dash import html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import os
import scanpy as sc
import json

from general import *

# must add this line in order for the app to be deployed successfully on Heroku
from app import server
from app import app
# import all pages in the app
from pages import page_home, page_analysis, page_h5ad, page_report

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
                    dbc.DropdownMenuItem("Home", href="/home"),                     
                    dbc.DropdownMenuItem("h5ad inspector", href="/h5ad"),                     
                    dbc.DropdownMenuItem("Graphical Analysis", href="/analysis"),                     
                    dbc.DropdownMenuItem("Report", href="/report"),                     
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

if os.path.exists(STARTPATH+"Raw_h5ad.sc"):
    load_analysis(STARTPATH+"Raw_h5ad.sc")
    pathname = "/analysis"
else:
    pathname = "/home"

# embedding the navigation bar
app.layout = html.Div([
    dcc.Location(id='url', pathname=pathname, refresh=False),
    navbar,
    html.Div(id='page-content')
])

@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):

    if pathname == '/home':

        return page_home.layout()

    elif pathname == '/h5ad':

        return page_h5ad.layout()
    
    elif pathname == '/analysis':

        return page_analysis.layout()

    elif pathname == '/report':

        return page_report.layout()


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