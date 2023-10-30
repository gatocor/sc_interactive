import dash
from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
import os
import scanpy as sc
import pandas as pd
from dash.exceptions import PreventUpdate
from shutil import rmtree

from general import *

from app import app

error_modal = dbc.Modal(
    [
        dbc.ModalHeader("Warning"),
        dbc.ModalBody(id="overlay-message",
                      children="An analysis from this file exist already."
                    ),
        dbc.ModalFooter([
            dbc.Button("Make new analysis", id="warning-new-button", className="ml-auto"),
            dbc.Button("Load existing analysis", id="warning-load-button", className="ml-auto")                         
                         ]),
    ],
    id="existing-analysis-modal",
    size="sm",
)

# change to app.layout if running as single page app instead
def layout():

    return dbc.Container([
        dcc.Dropdown(id="select-file",options=filesInFolderTree()),
        dbc.Button(id="load-button",children="Load"),
        html.Div(id="loaded-text", children=""),
        error_modal
    ])
    
def filesInFolderTree (startPath = STARTPATH, extension = ".h5ad", dirextension = '.sc'):

    listOfFileNames = [ os.path.join (root, name) \
                        for root, dirs, files in os.walk (startPath) \
                        for name in sorted(files) \
                        if name.endswith (extension) and not any( [dir.endswith(dirextension) for dir in root.split("/")] ) \
                      ] + \
                      [ os.path.join (root, name) \
                        for root, dirs, files in os.walk (startPath) \
                        for name in sorted(dirs) \
                        if name.endswith (dirextension) \
                      ]

    return listOfFileNames

@app.callback(
    dash.Output("loaded-text","children", allow_duplicate=True),
    dash.Output("select-file","options", allow_duplicate=True),
    dash.Output("existing-analysis-modal","is_open", allow_duplicate=True),
    dash.Input("load-button","n_clicks"),
    dash.State("select-file","value"),
    prevent_initial_call=True
)
def load(n_clicks, file):

    if n_clicks:

        if file:

            if file.endswith(ENDH5AD):

                folder = name_analysis_folder(file)

                if os.path.isdir(folder):

                    return "", filesInFolderTree(), True

                else:

                    create_analysis(folder, raw=file)
                    load_analysis(folder)

            elif file.endswith(ENDANALYSIS):

                load_analysis(file)

            return file, filesInFolderTree(), False

        else:

            raise PreventUpdate()
        
    raise PreventUpdate()

@app.callback(
    dash.Output("loaded-text","children", allow_duplicate=True),
    dash.Output("select-file","options", allow_duplicate=True),
    dash.Output("existing-analysis-modal","is_open", allow_duplicate=True),
    dash.Input("warning-new-button","n_clicks"),
    dash.State("select-file","value"),
    prevent_initial_call = True
)
def erase_and_start(n_clicks, file):

    if n_clicks:
        
        folder = name_analysis_folder(file)

        count = 1
        folder_new = f"{folder[:-len(ENDANALYSIS)]}_{count}{ENDANALYSIS}"
        while os.path.isdir(folder_new):
            count += 1
            folder_new = f"{folder}_{count}"

        create_analysis(folder_new, raw=file)
        load_analysis(folder_new)

        return file, filesInFolderTree(), False
    
    else:

        raise PreventUpdate()

@app.callback(
    dash.Output("loaded-text","children", allow_duplicate=True),
    dash.Output("select-file","options", allow_duplicate=True),
    dash.Output("existing-analysis-modal","is_open", allow_duplicate=True),
    dash.Input("warning-load-button","n_clicks"),
    dash.State("select-file","value"),
    prevent_initial_call = True
)
def erase_and_start(n_clicks, file):

    if n_clicks:
        
        folder = name_analysis_folder(file)

        load_analysis(folder)

        return file, filesInFolderTree(), False
    
    else:

        raise PreventUpdate()
