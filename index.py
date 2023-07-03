import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output, State
import dash_bootstrap_components as dbc
import os
import scanpy as sc

# must add this line in order for the app to be deployed successfully on Heroku
from app import server
from app import app
# import all pages in the app
from pages import config, home, quality_control

navbar = dbc.Navbar(
    dbc.Container(
        [
            html.A(
                # Use row and col to control vertical alignment of logo / brand
                dbc.Row(
                    [
                        dbc.Col(html.Img(src="/assets/scRNAseq.png", height="60px")),
                        dbc.Col(dbc.NavbarBrand("scRNAseq analysis", className="ml-2")),
                    ],
                    align="center",
                ),
                href="/home",
            ),
            dbc.Nav([
                dbc.DropdownMenu([
                    dbc.DropdownMenuItem("Quality Control", href="/quality_control")                        
                ],
                in_navbar = True,
                label = "Quality control",
                color="dark"
                ),
                dbc.DropdownMenu([
                    dbc.DropdownMenuItem("Quality Control", href="/quality_control")                        
                ],
                in_navbar = True,
                label = "Dimensionality Reduction",
                color="dark"
                ),
                dbc.DropdownMenu([
                    dbc.DropdownMenuItem("Clustering and Annotation", href="/quality_control")                        
                ],
                in_navbar = True,
                label = "Clustering",
                color="dark"
                ),
                dbc.DropdownMenu([
                    dbc.DropdownMenuItem("Clustering and Annotation", href="/quality_control")                        
                ],
                in_navbar = True,
                label = "Annotation",
                color="dark"
                ),
                dbc.DropdownMenu([
                    dbc.DropdownMenuItem("Pseudotime", href="/quality_control")                        
                ],
                in_navbar = True,
                label = "Time Integration",
                color="dark"
                )

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
    if pathname == '/quality_control':
        return quality_control.layout
    else:
        return home.layout

if __name__ == '__main__':
    app.run_server(debug=True)