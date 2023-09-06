import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from app import app

from .. import config

def pca_args(adata):

    options = []
    if "__interactive__" in adata.uns.keys():
        options = ["Raw"]+[i[2:] for i in adata.obsm.keys()] 
    
    return [
        "PCA",
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":None,
            "clearable":True,
            "options":options
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

def f_pca(adata, name_analysis, **kwargs):
        
    if kwargs["input"] == "Raw":
        adata_copy = sc.AnnData(X=adata.X)
    else:
        adata_copy = sc.AnnData(X=adata.obsm["X_"+kwargs["input"]])
    
    sc.pp.pca(adata_copy,
              n_comps=kwargs["n_comps"],
              zero_center=kwargs["zero_center"],
              svd_solver=kwargs["svd_solver"],
              random_state=kwargs["random_state"]
              )

    adata.obsm["X_"+name_analysis] = adata_copy.obsm["X_pca"]
    # adata.varm[name_analysis] = adata_copy.varm["PCs"]
    adata.uns[name_analysis] = adata_copy.uns["pca"]

    adata.uns["__interactive__"][name_analysis]["threshold"] = adata.obsm["X_"+name_analysis].shape[1]-1
    adata.uns["__interactive__"][name_analysis]["n_comps_plot"] = 3

def make_pca_plots1(adata, name_analysis):

    if "X_"+name_analysis in adata.obsm.keys():

        y = adata.uns[name_analysis]["variance_ratio"]
        threshold = adata.uns["__interactive__"][name_analysis]["threshold"]

        return [
            dbc.Col(
                dbc.Row([
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
                    dcc.Slider(id="pca_threshold",
                        min = 1,
                        max = len(y),
                        step = 1,
                        value = threshold 
                    ),
                ])
            ),
        ]
    else:
        return []

def make_pca_plots2(adata, name_analysis):

    if "X_"+name_analysis in adata.obsm.keys():

        n_comp = int(config.adata.uns["__interactive__"][name_analysis]["n_comps_plot"])

        fig = make_subplots(rows=n_comp-1, cols=n_comp-1, shared_yaxes=True, shared_xaxes=True)

        for i in range(n_comp-1):

            for j in range(i+1,n_comp):
                x_pca = adata.obsm["X_"+name_analysis][:,i]
                y_pca = adata.obsm["X_"+name_analysis][:,j]

                fig.add_trace(
                        go.Scattergl(
                                    x=x_pca,
                                    y=y_pca,
                                    mode='markers'    
                        ),
                        row=j, col=i+1
                )

        fig.update_layout(height=1200, width=1200, autosize=True, showlegend=False)
        # fig.update_xaxes(scaleanchor="y", scaleratio=1)
        # fig.update_yaxes(constrain='domain')

        return [
            dbc.Col(
                dbc.Row([
                        dcc.Slider(id="pca_n_comps_plot",
                                min = 2,
                                max = 7,
                                step = 1,
                                value = n_comp 
                            ),
                        dcc.Graph(
                            figure = fig
                        ),
                ]),
                width="auto"
            ),
        ]

    else:
        return []

@app.callback(
    dash.Output('pca_plots1', 'children', allow_duplicate=True),
    dash.Input('pca_threshold', 'value'),
    dash.State('dimensionality_reduction_method','children'),
    prevent_initial_call=True
)
def pca_threshold(pca_threshold, name_analysis):

    config.adata.uns["__interactive__"][name_analysis]["threshold"] = pca_threshold

    return  make_pca_plots1(config.adata, name_analysis)

@app.callback(
    dash.Output('pca_plots2', 'children', allow_duplicate=True),
    dash.Input('pca_n_comps_plot', 'value'),
    dash.State('dimensionality_reduction_method','children'),
    prevent_initial_call=True
)
def pca_n_comps(pca_n_comps, name_analysis):

    config.adata.uns["__interactive__"][name_analysis]["n_comps_plot"] = pca_n_comps

    return make_pca_plots2(config.adata, name_analysis)