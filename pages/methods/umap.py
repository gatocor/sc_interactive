import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash
import scrublet
from scipy.stats import mode

from ..functions import *
from ..graph import *
from ..plots import *

from app import app

from .. import config

def args_umap():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_batch = get_batch_keys()

    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Use the indicated neighbors.",
            "value":None,
            "clearable":False,
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
            "input":"Input",
            "name":"min_dist",
            "description":"The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out. The default of in the umap-learn package is 0.1.",
            "value":.5,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"spread",
            "description":"The effective scale of embedded points. In combination with min_dist this determines how clustered/clumped the embedded points are.",
            "value":1.0,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"n_dimensions",
            "description":"The number of dimensions of the embedding.",
            "value":2,
            "clearable":False,
            "options":[2,3]
        },
        {
            "input":"Input",
            "name":"maxiter",
            "description":"The number of iterations (epochs) of the optimization. Called n_epochs in the original UMAP.",
            "value":None,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"alpha",
            "description":"The initial learning rate for the embedding optimization.",
            "value":1.0,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"gamma",
            "description":"Weighting applied to negative samples in low dimensional embedding optimization. Values higher than one will result in greater weight being given to negative samples.",
            "value":1.0,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"negative_sample_rate",
            "description":"The number of negative edge/1-simplex samples to use per positive edge/1-simplex sample in optimizing the low dimensional embedding.",
            "value":5,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"init_pos",
            "description":"How to initialize the low dimensional embedding.",
            "value":'spectral',
            "clearable":False,
            "options":['spectral','random']
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"If int, random_state is the seed used by the random number generator; If RandomState or Generator, random_state is the random number generator; If None, the random number generator is the RandomState instance used by np.random.",
            "value":0,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"a",
            "description":"More specific parameters controlling the embedding. If None these values are set automatically as determined by min_dist and spread.",
            "value":None,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"b",
            "description":"More specific parameters controlling the embedding. If None these values are set automatically as determined by min_dist and spread.",
            "value":None,
            "type":"number"
        },
    ]

def f_umap(name_analysis, kwargs, sub_name, sub):

    node = get_node(kwargs["input"])
    x = config.adata.uns[kwargs["input"]+"_"+str(sub_name)]["connectivities"]
    adata_copy = sc.AnnData(X = np.zeros([x.shape[0],2]))
    adata_copy.uns["neighbors"] = {"connectivities_key":"connectivities","params":{"method":node["data"]["parameters"]["method"]}}
    adata_copy.obsp["connectivities"] = x

    sc.tl.umap(adata_copy,
                neighbors_key=kwargs["input"],
                min_dist=kwargs["min_dist"],
                spread=kwargs["spread"],
                n_components=kwargs["n_dimensions"],
                maxiter=kwargs["maxiter"],
                alpha=kwargs["alpha"],
                gamma=kwargs["gamma"],
                negative_sample_rate=kwargs["negative_sample_rate"],
                init_pos=kwargs["init_pos"],
                random_state=kwargs["random_state"],
                a=kwargs["a"],
                b=kwargs["b"],
                )

    pos = get_node_pos(name_analysis)
    if "obsm" not in config.graph[pos]["data"]:
        config.graph[pos]["data"]["obsm"] = np.zeros([config.adata.shape[0],kwargs["n_dimensions"]])
    elif config.graph[pos]["data"]["obsm"].shape[1] != kwargs["n_dimensions"]:
        config.graph[pos]["data"]["obsm"] = np.zeros([config.adata.shape[0],kwargs["n_dimensions"]])

    config.graph[pos]["data"]["obsm"][sub,:] = adata_copy.obsm["X_umap"]
    config.graph[pos]["data"]["plot"] = {"color":None}

def rm_umap(name_analysis):

    return

def rename_umap(name_analysis, name_new_analysis):

    return

def plot_umap(name_analysis):

    node = get_node(name_analysis)
    if not node["data"]["computed"]:
        return []    

    X = np.array(node["data"]["obsm"])
    c = get_observable(node["data"]["plot"]["color"])
    if type(c) == type(None):
        c = np.ones(X.shape[0],dtype=bool)

    bs = np.array(get_observable(node["data"]["parameters"]["batch_key"]))
    if bs == None:
        bs = np.ones(X.shape[0],dtype=bool)

    l = [
        dcc.Dropdown(id="umap_color",
                value = node["data"]["plot"]["color"],
                options = list_observables()
            ),
    ]
    if node["data"]["parameters"]["n_dimensions"] == 2:

        for b in np.unique(bs):

            sub = bs == b

            x = X[sub,0]
            y = X[sub,1]

            l += [
                dbc.Col(),
                dbc.Col(
                    dcc.Graph(
                        figure = {'data':[
                                    go.Scattergl(
                                        x=x,
                                        y=y,
                                        mode='markers',
                                        marker={"color":qualitative_colors(c)}
                                    )],
                                    'layout':{
                                            'yaxis': {
                                                'scaleanchor': 'x',
                                                'scaleratio': 1
                                            },
                                            'width':900,
                                            'height':800,
                                    },
                        }
                    )
                ),
                dbc.Col()
            ]
    else:

        for b in np.unique(bs):

            sub = bs == b

            x = X[sub,0]
            y = X[sub,1]
            z = X[sub,2]

            l += [
                dbc.Col(),
                dbc.Col(
                    dcc.Graph(
                        figure = {'data':[
                                    go.Scatter3d(
                                        x=x,
                                        y=y,
                                        z=z,
                                        mode='markers',
                                        marker={"color":qualitative_colors(c)},
                                    )],
                                'layout':{
                                        'yaxis': {
                                            'scaleanchor': 'x',
                                            'scaleratio': 1
                                        },
                                        'zaxis': {
                                            'scaleanchor': 'x',
                                            'scaleratio': 1
                                        },
                                        'width':900,
                                        'height':800,
                                },
                        }
                    )
                ),
                dbc.Col()
            ]
    
    return l

@app.callback(
    dash.Output("analysis_plot", "children", allow_duplicate=True),
    dash.Input("umap_color", "value"),
    prevent_initial_call=True
)
def umap_threshold(pca_color):

    prevent_race("umap")

    pos = get_node_pos(config.selected)
    config.graph[pos]["data"]["plot"]["color"] = pca_color

    return  plot_umap(config.selected)