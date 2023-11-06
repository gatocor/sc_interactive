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

from general import *

tsne_args = {

    "execution" : [
        ARGINPUT,

        {
            "input":"Dropdown",
            "name":"use_rep",
            "description":"Use the indicated representation. 'X' or any key for .obsm is valid. If None, the representation is chosen automatically: For .n_vars < scanpy.settings.N_PCS (default: 50), .X is used, otherwise ‘X_pca’ is used. If ‘X_pca’ is not present, it’s computed with default parameters or n_pcs if present.",
            "properties" : {
                "value":"X",
                "options":{"function":"['X']+[i for i in config.adata.obsm.keys()]"}
            }
        },
        {
            "input":"Dropdown",
            "name":"n_pcs",
            "description":"Use this many components.",
            "properties" : {
                "value":None,
                "options":{"function":"neighbors_components()"}
            },
            "recomputeAfter":["use_rep"]
        },
        {
            "input":"Input",
            "name":"perplexity",
            "description":"The perplexity is related to the number of nearest neighbors that is used in other manifold learning algorithms. Larger datasets usually require a larger perplexity. Consider selecting a value between 5 and 50. The choice is not extremely critical since t-SNE is quite insensitive to this parameter.",
            "properties" : {
                "value":30,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"metric",
            "description":"Distance metric calculate neighbors on.",
            "properties" : {
                "value":"euclidean",
                "options":['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan','braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
            }
        },
        {
            "input":"Input",
            "name":"early_exaggeration",
            "description":"Controls how tight natural clusters in the original space are in the embedded space and how much space will be between them. For larger values, the space between natural clusters will be larger in the embedded space. Again, the choice of this parameter is not very critical. If the cost function increases during initial optimization, the early exaggeration factor or the learning rate might be too high.",
            "properties" : {
                "value":12,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"learning_rate",
            "description":"Note that the R-package “Rtsne” uses a default of 200. The learning rate can be a critical parameter. It should be between 100 and 1000. If the cost function increases during initial optimization, the early exaggeration factor or the learning rate might be too high. If the cost function gets stuck in a bad local minimum increasing the learning rate helps sometimes.",
            "properties" : {
                "value":1000,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change this to use different intial states for the optimization. If None, the initial state is not reproducible.",
            "properties" : {
                "value":0,
                "type":"number"
            }
        },
    ],

    "postexecution" : [],

    "plot" : ARGS_COLOR

}

def tsne_f(adata, kwargs):

    sc.tl.tsne(adata,
        n_pcs = kwargs["n_pcs"],
        use_rep = kwargs["use_rep"],
        perplexity = kwargs["perplexity"],
        early_exaggeration = kwargs["early_exaggeration"],
        learning_rate = kwargs["learning_rate"],
    )

def tsne_plot():

    c = get_color()
    X = config.adata.obsm["X_tsne"]

    fig = px.scatter(
                    x=X[:,0],
                    y=X[:,1],
                    color=c,
                    height=PLOTHEIGHT,
                    width=PLOTWIDTH
            )
    
    return  plot_center(dcc.Graph(figure=fig))

config.methods["tsne"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": tsne_args,

    "function": tsne_f,

    "plot": tsne_plot,

}