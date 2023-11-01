import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
from dash import dcc, dash_table
import plotly.graph_objs as go
from plotly.subplots import make_subplots

from general import *

args_pca = {

    "execution" : [
        ARGINPUT,
        {
            "input":"BooleanSwitch",
            "name":"use_highly_varying",
            "description":"Whether to use highly variable genes only.",
            "properties":{
                "on":True,
            }
        },
        {
            "input":"Input",
            "name":"n_comps",
            "description":"Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation.",
            "properties":{
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"BooleanSwitch",
            "name":"zero_center",
            "description":"If True, compute standard PCA from covariance matrix. If False, omit zero-centering variables (uses TruncatedSVD), which allows to handle sparse input efficiently. Passing None decides automatically based on sparseness of the data.",
            "properties":{
                "value":True,
            }
        },
        {
            "input":"Dropdown",
            "name":"svd_solver",
            "description":"Efficient computation of the principal components of a sparse matrix currently only works with the 'arpack' or 'lobpcg' solvers.",
            "properties":{
                "value":"arpack",
                "clearable":False,
                "options":["arpack", "randomized", "auto","lobpcg"] 
            }
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change to use different initial states for the optimization.",
            "properties":{
                "value":0,
                "type":"number"
            }
        },
    ],

    "postexecution" : [],

    "plot" : [

    ]
}

def pca_f(adata, kwargs):
        
    sc.pp.pca(config.adata,
              use_highly_varying=kwargs["use_highly_varying"],
              n_comps=kwargs["n_comps"],
              zero_center=kwargs["zero_center"],
              svd_solver=kwargs["svd_solver"],
              random_state=kwargs["random_state"]
              )

def pca_plot():

    if not node["data"]["computed"]:
        return []

    columns = node["data"]["threshold"]["columns"]
    data = node["data"]["threshold"]["data"]
    color = node["data"]["plot"]["color"]

    l = [
            dbc.Row([
                dash_table.DataTable(id="pca_threshold",
                    columns= columns,
                    data = data 
                ),
            ]),
            dcc.Dropdown(id="pca_color",
                    value = color,
                    options = list_observables()
                ),
    ]

    for b,c_sim in config.graph[pos]["data"]["uns"].items():

        if b == "null":
            b = None

        if b == None:
            sub = np.ones(np.array(config.graph[pos]["data"]["obsm"]["PCA"]).shape[0],dtype=bool)
        else:
            sub = np.array(config.graph[pos]["data"]["batch"]) == b

        y = c_sim["variance_ratio"]

        threshold = int(get_table_value(data,b,"pca maxPCA"))
        n_comp = int(get_table_value(data,b,"pca nPlots"))
        c = get_observable(color)

        fig = make_subplots(rows=n_comp-1, cols=n_comp-1, shared_yaxes=True, shared_xaxes=True)

        X = np.array(config.graph[pos]["data"]["obsm"]["PCA"])

        for i in range(n_comp-1):

            for j in range(i+1,n_comp):
                x_pca = X[sub,i]
                y_pca = X[sub,j]

                fig.add_trace(
                        go.Scattergl(
                                    x=x_pca,
                                    y=y_pca,
                                    mode="markers",
                                    marker=dict(color=qualitative_colors(c)),    
                        ),
                        row=j, col=i+1
                )

        fig.update_layout(height=1200, width=1200, autosize=True, showlegend=False)

        l += [
            dbc.Row([
                dcc.Graph(
                    figure = {"data":[
                                    go.Scatter(
                                        x=np.array(range(len(y))),
                                        y=y,
                                        mode="markers",
                                        name="Variance ratio",
                                    ),
                                    go.Scatter(
                                        x=[threshold+.5, threshold+.5],
                                        y=[0,max(y)],
                                        mode="lines",
                                        name="Threshold",
                                        line=dict(color="red")
                                    )
                                ]
                            }
                ),
            ]),
            dbc.Row([
                dbc.Col(),
                dbc.Col(
                    dcc.Graph(
                        figure = fig
                    ),
                ),
                dbc.Col(),
            ]),
        ]

    return l