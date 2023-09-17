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
    options_higly_variable = [i["data"]["name"] for i in get_nodes() if i["data"]["method"] in ["feature_selection"]]
    options_batch = get_batch_keys()
    
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
            "name":"batch_key",
            "description":"str, optional (default: None) Batch key to use. The highly varying will be computed independently in each set of cells separated by batch. If None, use the full dataset.",
            "value":None,
            "clearable":True,
            "options":options_batch,
            "summary":True
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
            "description":"Efficient computation of the principal components of a sparse matrix currently only works with the 'arpack' or 'lobpcg' solvers.",
            "value":"arpack",
            "clearable":False,
            "options":["arpack", "randomized", "auto","lobpcg"] 
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change to use different initial states for the optimization.",
            "value":0,
            "type":"number"
        },
    ]

def f_pca(name_analysis, kwargs, sub_name, sub):
        
    node = get_node(name_analysis)
    pos = get_node_pos(name_analysis)

    adata_copy = sc.AnnData(X=config.adata.X[sub,:])
    if kwargs["use_highly_varying"] != None:
        node_input = get_node(kwargs["use_highly_varying"])
        adata_copy.var["highly_variable"] = node_input["data"]["var"]["highly_variable_"+str(sub_name)]

    if kwargs["n_comps"] == None:
        kwargs["n_comps"] = 50

    if "obsm" not in node["data"].keys():
        config.graph[pos]["data"]["obsm"] = {"PCA":np.zeros([config.adata.shape[0],kwargs["n_comps"]])}
        config.graph[pos]["data"]["uns"] = {}

    sc.pp.pca(adata_copy,
              n_comps=kwargs["n_comps"],
              zero_center=kwargs["zero_center"],
              svd_solver=kwargs["svd_solver"],
              random_state=kwargs["random_state"]
              )

    config.graph[pos]["data"]["obsm"]["PCA"][sub,:] = adata_copy.obsm["X_pca"]
    config.graph[pos]["data"]["uns"][sub_name] =  {
        "variance":adata_copy.uns["pca"]["variance"],
        "variance_ratio":adata_copy.uns["pca"]["variance"],
    }

    config.graph[pos]["data"]["plot"] = {
        "color": None
    }

    add = ["maxPCA","nPlots"]
    if kwargs["batch_key"] == None:
        rows = [None]
    else:
        rows = np.sort(config.adata.obs[kwargs["batch_key"]].unique())
    
    columns, data = make_thresholds_table(["pca"], rows, add)

    #Fill table 
    for i in rows:
        set_table_value(data, i, "pca maxPCA", kwargs["n_comps"])
        set_table_value(data, i, "pca nPlots", 3)

    config.graph[pos]["data"]["threshold"] = {"columns":columns,"data":data}


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
    pos = get_node_pos(name_analysis)

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

@app.callback(
    dash.Output("analysis_plot", "children", allow_duplicate=True),
    dash.Input("pca_threshold", "data"),
    dash.Input("pca_color", "value"),
    prevent_initial_call=True
)
def pca_threshold(pca_threshold, pca_color):

    prevent_race("pca")

    pos = get_node_pos(config.selected)
    config.graph[pos]["data"]["threshold"]["data"] = pca_threshold
    config.graph[pos]["data"]["plot"]["color"] = pca_color

    return  plot_pca(config.selected)

@app.callback(
    dash.Output("analysis_batch_key", "value", allow_duplicate=True),
    dash.Input("analysis_use_highly_varying", "value"),
    dash.Input("analysis_batch_key", "value"),
    prevent_initial_call=True
)
def match_highly_varying(hvgs, bk):

    prevent_race("pca")

    if hvgs == None:
        return bk
    else:
        return get_node(hvgs)["data"]["parameters"]["batch_key"]
