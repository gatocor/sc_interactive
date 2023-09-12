import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
from dash import dcc, dash_table
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from ..functions import *

from app import app

from .. import config

def args_pca():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_higly_variable = [i for i in config.adata.var.columns if config.adata.var.dtypes[i] in [bool]]
    
    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":None,
            "clearable":True,
            "options":options
        },
        {
            "input":"Dropdown",
            "name":"use_highly_varying",
            "description":"Whether to use highly variable genes only.",
            "value":None,
            "clearable":True,
            "options":options_higly_variable
        },
        {
            "input":"Input",
            "name":"n_comps",
            "description":"Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation.",
            "value":None,
            "type":"number"
        },
        {
            "input":"BooleanSwitch",
            "name":"zero_center",
            "description":"If True, compute standard PCA from covariance matrix. If False, omit zero-centering variables (uses TruncatedSVD), which allows to handle sparse input efficiently. Passing None decides automatically based on sparseness of the data.",
            "value":True,
        },
        {
            "input":"Dropdown",
            "name":"svd_solver",
            "description":"Efficient computation of the principal components of a sparse matrix currently only works with the 'arpack’ or 'lobpcg' solvers.",
            "value":'arpack',
            "clearable":False,
            "options":['arpack', 'randomized', 'auto','lobpcg'] 
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change to use different initial states for the optimization.",
            "value":0,
            "type":"number"
        },
    ]

def f_pca(name_analysis, kwargs):
        
    adata_copy = sc.AnnData(X=config.adata.X)
    if kwargs['use_highly_varying'] != None:
        adata_copy.var["higly_variable"] = config.adata.var[kwargs['use_highly_varying']].values
    
    sc.pp.pca(adata_copy,
              n_comps=kwargs["n_comps"],
              zero_center=kwargs["zero_center"],
              svd_solver=kwargs["svd_solver"],
              random_state=kwargs["random_state"]
              )

    config.adata.obsm[name_analysis] = adata_copy.obsm["X_pca"]
    # adata.varm[name_analysis] = adata_copy.varm["PCs"]
    config.adata.uns[name_analysis] = adata_copy.uns["pca"]

    pos = get_node_pos(name_analysis)
    config.graph[pos]['data']["threshold"] = config.adata.obsm[name_analysis].shape[1]-1
    config.graph[pos]['data']["n_comps_plot"] = 3
    config.graph[pos]['data']["color"] = None

def rm_pca(name_analysis):

    del config.adata.obsm[name_analysis]
    del config.adata.uns[name_analysis]

def rename_pca(name_analysis, name_new_analysis):

    config.adata.obsm[name_new_analysis] = config.adata.obsm[name_analysis]
    config.adata.uns[name_new_analysis] = config.adata.uns[name_analysis]
    del config.adata.obsm[name_analysis]
    del config.adata.uns[name_analysis]

def plot_pca(name_analysis):

    node = get_node(name_analysis)
    if not node['data']['computed']:
        return []

    y = config.adata.uns[name_analysis]["variance_ratio"]

    threshold = node['data']["threshold"]
    n_comp = int(node['data']["n_comps_plot"])
    color = node['data']["color"]
    if color == None:
        c = np.zeros_like(config.adata.obsm[name_analysis][:,0])
    else:
        c = config.adata.obs[color].values

    fig = make_subplots(rows=n_comp-1, cols=n_comp-1, shared_yaxes=True, shared_xaxes=True)

    for i in range(n_comp-1):

        for j in range(i+1,n_comp):
            x_pca = config.adata.obsm[name_analysis][:,i]
            y_pca = config.adata.obsm[name_analysis][:,j]

            fig.add_trace(
                    go.Scattergl(
                                x=x_pca,
                                y=y_pca,
                                mode='markers',
                                marker=dict(color=qualitative_colors(c)),    
                    ),
                    row=j, col=i+1
            )

    fig.update_layout(height=1200, width=1200, autosize=True, showlegend=False)

    return [
        dbc.Row([
            dcc.Slider(id="pca_threshold",
                min = 1,
                max = len(y),
                step = 1,
                value = threshold 
            ),
            dcc.Graph(
                figure = {'data':[
                                go.Scatter(
                                    x=np.array(range(len(y))),
                                    y=y,
                                    mode='markers',
                                    name='Variance ratio',
                                ),
                                go.Scatter(
                                    x=[threshold+.5, threshold+.5],
                                    y=[0,max(y)],
                                    mode='lines',
                                    name='Threshold',
                                    line=dict(color='red')
                                )
                            ]
                        }
            ),
        ]),
        dbc.Row([
                dcc.Slider(id="pca_n_comps_plot",
                        min = 2,
                        max = 7,
                        step = 1,
                        value = n_comp 
                    ),
                dcc.Dropdown(id="pca_color",
                        value = color,
                        options = config.adata.obs.columns
                    ),
                dcc.Graph(
                    figure = fig
                ),
        ]),
    ]

@app.callback(
    dash.Output('analysis_plot', 'children', allow_duplicate=True),
    dash.Input('pca_threshold', 'value'),
    dash.Input('pca_n_comps_plot', 'value'),
    dash.Input('pca_color', 'value'),
    prevent_initial_call=True
)
def pca_threshold(pca_threshold, pca_n_comps, pca_color):

    prevent_race("pca")

    pos = get_node_pos(config.selected)
    config.graph[pos]['data']["threshold"] = pca_threshold
    config.graph[pos]['data']["n_comps_plot"] = pca_n_comps
    config.graph[pos]['data']["color"] = pca_color

    return  plot_pca(config.selected)