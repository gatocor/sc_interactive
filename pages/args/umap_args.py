import numpy as np

def umap_args(adata):
    return [
        "UMAP",
        {
            "input":"BooleanSwitch",
            "name":"log1p",
            "description":"If True, log scale the count matrix before PCA analysis.",
            "value":True,
        },
        {
            "input":"Dropdown",
            "name":"highly_variable",
            "description":"Genes to use to for the dimensionality reduction.",
            "value":None,
            "clearable":True,
            "options":[i for i,j in adata.uns["dimensionality_reduction"].items() if j["type"]=="Feature Selection"] 
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
