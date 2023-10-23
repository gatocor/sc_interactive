import dash
from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
import os
import scanpy as sc
from . import config
import pandas as pd
from .functions import *
import dash_editor_components
import re
from .graph import save_report

from app import app

from PIL import Image

load_analysis("../Raw_h5ad.sc")

# change to app.layout if running as single page app instead
def layout():

    return dbc.Container([
        dbc.Row([
            dbc.Col(
                dbc.Switch(id="show-editor",value=False,label="Activate editor")
            ),
            dbc.Col(
                dbc.Button(id="save-report",children="Save",disabled=False)
            ),
        ]),
        dbc.Row(
            id = "report",
            children = [
                dbc.Col(id="markdown"),
                dbc.Col(
                    dash_editor_components.PythonEditor(
                        id='editor',
                        value=config.report
                    )
                ),
            ]
        ),
        html.Div(id="dump",children="")
    ])

@app.callback(
    dash.Output("report","children"),
    dash.Output("save-report","disabled"),
    dash.Input("show-editor","value"),
    dash.State("report","children"),
)
def activate_report(val, s):
        
    if val != None:

        l = markdown_to_dash(config.report)

        if val:

            d = [
                    dbc.Col(id="markdown", children=l),
                    dbc.Col(
                        dash_editor_components.PythonEditor(
                            id='editor',
                            value=config.report
                        )
                    ),
                ]
            
            active = False

        else:

            d = [
                    dbc.Col(id="markdown", children=l),
                ]

            active = True

        return d, active

    else:

        raise PreventUpdate() 

@app.callback(
    dash.Output("dump","children"),
    dash.Input("save-report","n_clicks")
)
def editor(value):

    if value:

        save_report()

        return ""
    
    else:

        raise PreventUpdate()

@app.callback(
    dash.Output("markdown","children"),
    dash.Input("editor","value")
)
def editor(value):

    if value:
        config.report = value

        l = markdown_to_dash(value)

        return l
    
    else:

        raise PreventUpdate()

def markdown_to_dash(value):
    
    res = re.findall(r'\!\[.*\]\(.*?\)', value)

    l = []
    for i in res:
        pre, value = value.split(i)
        l += [dcc.Markdown(children=pre)]
        img_source = i.split("](")[1][:-1]
        img = Image.open(config.analysis_folder+"/report/"+img_source)
        l += [html.Img(
                src=img,
                width="100%")]
        
    l += [dcc.Markdown(children=value)]
        
    return l