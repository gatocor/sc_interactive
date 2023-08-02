import numpy as np

def neighbors_args(adata):

    options = []
    if "dimensionality_reduction" in adata.uns.keys():
        options = [i for i,j in adata.uns["dimensionality_reduction"].items() if j["type"]=="Feature Selection"] 

    return [
        "Neihgbors",
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Use the indicated representation. 'X' or any key for .obsm is valid. If None, the representation is chosen automatically: For .n_vars < 50, .X is used, otherwise 'X_pca' is used. If 'X_pca' is not present, it's computed with default parameters.",
            "value":None,
            "clearable":True,
            "options":options
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
