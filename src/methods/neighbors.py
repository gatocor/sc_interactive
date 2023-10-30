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

def args_neighbors():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_batch = get_batch_keys()

    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Use the indicated representation. 'X' or any key for .obsm is valid. If None, the representation is chosen automatically: For .n_vars < 50, .X is used, otherwise 'X_pca' is used. If 'X_pca' is not present, it's computed with default parameters.",
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
            "input":"Input",
            "name":"n_components",
            "description":"Use that many graphs from the input representation.",
            "value":15,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"n_neighbors",
            "description":"The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100. If knn is True, number of nearest neighbors to be searched. If knn is False, a Gaussian kernel width is set to the distance of the n_neighbors neighbor.",
            "value":15,
            "type":"number"
        },
        {
            "input":"BooleanSwitch",
            "name":"knn",
            "description":"If True, use a hard threshold to restrict the number of neighbors to n_neighbors, that is, consider a knn graph. Otherwise, use a Gaussian Kernel to assign low weights to neighbors more distant than the n_neighbors nearest neighbor.",
            "value":True,
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"A numpy random seed.",
            "value":0,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"method",
            "description":"Use 'umap' [McInnes18] or 'gauss' (Gauss kernel following [Coifman05] with adaptive width [Haghverdi16]) for computing connectivities. Use 'rapids' for the RAPIDS implementation of UMAP (experimental, GPU only).",
            "value":'umap',
            "clearable":False,
            "options":['umap', 'gauss', 'rapids']
        },
        {
            "input":"Dropdown",
            "name":"metric",
            "description":"A known metric's name or a callable that returns a distance.",
            "value":'correlation',
            "clearable":False,
            "options":['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan','braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
        },
    ]

def f_neighbors(name_analysis, kwargs, sub_name, sub):
        
    node_input = get_node(kwargs["input"])
    if "obsm" in node_input["data"].keys():
        n = list(node_input["data"]["obsm"])
        adata_copy = sc.AnnData(X=np.array(node_input["data"]["obsm"][n[0]])[sub,:kwargs["n_components"]])
    else:
        adata_copy = sc.AnnData(X=config.adata.X)

    sc.pp.neighbors(adata_copy,
              n_neighbors=kwargs["n_neighbors"],
              knn=kwargs["knn"],
              random_state=kwargs["random_state"],
              method=kwargs["method"],
              metric=kwargs["metric"],
    )

    if name_analysis not in config.adata.uns.keys():
        config.adata.uns[name_analysis] = {}

    config.adata.uns[name_analysis+"_"+str(sub_name)] = {
            "connectivities":adata_copy.obsp["connectivities"],
            "distances":adata_copy.obsp["distances"]
        }

def rm_neighbors(name_analysis):

    del config.adata.uns[name_analysis]
    return

def rename_neighbors(name_analysis, name_new_analysis):

    config.adata.uns[name_new_analysis] = config.adata.uns[name_analysis]
    del config.adata.uns[name_analysis]

    return

def plot_neighbors(name_analysis):

    return []