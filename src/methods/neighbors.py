import scanpy as sc
from general import *

neighbors_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"Input",
            "name":"n_neighbors",
            "description":"The size of local neighborhood (in terms of number of neighboring data points) used for manifold approximation. Larger values result in more global views of the manifold, while smaller values result in more local data being preserved. In general values should be in the range 2 to 100. If knn is True, number of nearest neighbors to be searched. If knn is False, a Gaussian kernel width is set to the distance of the n_neighbors neighbor.",
            "properties" : {
                "value":15,
                "type":"number"
            }
        },
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
            "input":"BooleanSwitch",
            "name":"knn",
            "description":"If True, use a hard threshold to restrict the number of neighbors to n_neighbors, that is, consider a knn graph. Otherwise, use a Gaussian Kernel to assign low weights to neighbors more distant than the n_neighbors nearest neighbor.",
            "properties" : {
                "on":True,
            }
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"A numpy random seed.",
            "properties" : {
                "value":0,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"method",
            "description":"Use 'umap' [McInnes18] or 'gauss' (Gauss kernel following [Coifman05] with adaptive width [Haghverdi16]) for computing connectivities. Use 'rapids' for the RAPIDS implementation of UMAP (experimental, GPU only).",
            "properties" : {
                "value":'umap',
                "clearable":False,
                "options":['umap', 'gauss', 'rapids']
            }
        },
        {
            "input":"Dropdown",
            "name":"metric",
            "description":"A known metric's name or a callable that returns a distance.",
            "properties" : {
                "value":'correlation',
                "clearable":False,
                "options":['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan','braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule']
            }
        },
    ],

    "postexecution" : [],

    "plot" : []

}

def neighbors_components():
    if config.active_node_parameters["use_rep"] == "X":
        return np.arange(1,config.adata.X.shape[1])
    else:
        return np.arange(1,config.adata.obsm[config.active_node_parameters["use_rep"]].shape[1])

def neighbors_f(adata, kwargs):
        
    sc.pp.neighbors(adata,
              n_neighbors=kwargs["n_neighbors"],
              knn=kwargs["knn"],
              random_state=kwargs["random_state"],
              method=kwargs["method"],
              metric=kwargs["metric"],
    )

def neighbors_plot():

    return []

config.methods["neighbors"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": neighbors_args,

    "function": neighbors_f,

    "plot": neighbors_plot,

}