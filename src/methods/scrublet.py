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

args = {

    "execution":[
        ARGINPUT,
        ARGBATCH,
        {
            "input":"BooleanSwitch",
            "name":"qc_before_computation",
            "description":"bool, optional (default: True) Remove cells from other quality control measures before computing the metrics.",
            "value":True,
        },
        {
            "input":"Input",
            "name":"synthetic_doublet_umi_subsampling",
            "description":"float, optional (default: 1.0) Rate for sampling UMIs when creating synthetic doublets. If 1.0, each doublet is created by simply adding the UMIs from two randomly sampled observed transcriptomes. For values less than 1, the UMI counts are added and then randomly sampled at the specified rate.",
            "value":1.0,
            "type":"number"
        },
        {
            "input":"BooleanSwitch",
            "name":"use_approx_neighbors",
            "description":"bool, optional (default: True) Use approximate nearest neighbor method (annoy) for the KNN classifier.",
            "value":True,
        },
        {
            "input":"Dropdown",
            "name":"distance_metric",
            "description":"str, optional (default: 'euclidean') Distance metric used when finding nearest neighbors. For list of valid values, see the documentation for annoy (if 'use_approx_neighbors' is True) or sklearn.neighbors.NearestNeighbors (if 'use_approx_neighbors' is False).",
            "value":"euclidean",
            "clearable":False,
            "options":["euclidean","manhattan","angular","hamming","dot"],
            "summary":True
        },
        {
            "input":"Input",
            "name":"min_counts",
            "description":"float, optional (default: 3) Used for gene filtering prior to PCA. Genes expressed at fewer than  'min_counts' in fewer than 'min_cells' (see below) are excluded.",
            "value":3,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"min_cells",
            "description":"int, optional (default: 3) Used for gene filtering prior to PCA. Genes expressed at fewer than  'min_counts' (see above) in fewer than 'min_cells' are excluded.",
            "value":3,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"min_gene_variability_pctl",
            "description":"float, optional (default: 85.0) Used for gene filtering prior to PCA. Keep the most highly variable genes (in the top min_gene_variability_pctl percentile), as measured by  the v-statistic [Klein et al., Cell 2015].",
            "value":85.0,
            "type":"number"
        },
        {
            "input":"BooleanSwitch",
            "name":"log_transform",
            "description":"bool, optional (default: False) If True, log-transform the counts matrix (log10(1+TPM)).  'sklearn.decomposition.TruncatedSVD' will be used for dimensionality reduction, unless 'mean_center' is True.",
            "value":False,
            "summary":True
        },
        {
            "input":"BooleanSwitch",
            "name":"mean_center",
            "description":"bool, optional (default: True) If True, center the data such that each gene has a mean of 0. 'sklearn.decomposition.PCA' will be used for dimensionality reduction.",
            "value":True,
        },
        {
            "input":"BooleanSwitch",
            "name":"normalize_variance",
            "description":"bool, optional (default: True) If True, normalize the data such that each gene has a variance of 1. 'sklearn.decomposition.TruncatedSVD' will be used for dimensionality reduction, unless 'mean_center' is True.",
            "value":True,
        },
        {
            "input":"Input",
            "name":"n_prin_comps",
            "description":"int, optional (default: 30) Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction.",
            "value":30,
            "type":"number",
            "summary":True
        },
        {
            "input":"Dropdown",
            "name":"svd_solver",
            "description":"str, optional (default: 'arpack') SVD solver to use. See available options for  'svd_solver' from 'sklearn.decomposition.PCA' or 'algorithm' from 'sklearn.decomposition.TruncatedSVD'",
            "value": "arpack",
            "clearable":False,
            "options": ["arpack"]
        }
    ],

    "postexecution" : [
        {
            "input":"AgTable",
            "name":"thresholds",
            "description":"Thresholds to apply to the data",
            "header":[
                { "headerName": "Batch", "field":"batch", "editable": False },
                { "headerName": "Max", "field":"max", "editable": True },
            ],
            "value":{"function":"scrublet_data()"},
        },
    ],

    "plot" : [
        {
            "input":"BooleanSwitch",
            "name":"show_scores",
            "description":"bool, optional (default: True) Pot scores of imputed cells.",
            "value":True,
        },
    ]

}

def scrublet_data():

    selected = config.selected
    adata = config.adata

    parameters = adata.uns[selected]["parameters"]

    batch = [""]
    if parameters["batch"] != None:
        batch = np.unique(config.adata.obs[parameters["batch"]].values)

    data = []
    j = f"{config.selected}--doublet_scores_obs"
    for b in batch:
        data.append({"batch":str(b),"max":str(adata.obs[j].max())})

    return data

def scrublet_(X, kwargs):

    scrub = scrublet.Scrublet(X)

    scrub.scrub_doublets(
        synthetic_doublet_umi_subsampling = kwargs["synthetic_doublet_umi_subsampling"], 
        use_approx_neighbors = kwargs["use_approx_neighbors"], 
        distance_metric = kwargs["distance_metric"], 
        min_counts = kwargs['min_counts'],
        min_cells = kwargs["min_cells"], 
        min_gene_variability_pctl = kwargs["min_gene_variability_pctl"], 
        log_transform = kwargs["log_transform"], 
        mean_center = kwargs["mean_center"], 
        normalize_variance = kwargs["normalize_variance"], 
        n_prin_comps = kwargs["n_prin_comps"], 
        svd_solver= kwargs["svd_solver"]
    )
            
    X = scrublet.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)

    return X, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_

def scrublet_f(adata, inputArgs, kwargs):

    pos = get_node_pos(config.selected)

    X = adata.X

    X_umap = np.zeros((X.shape[0],2))
    doublet_scores_obs = np.zeros(X.shape[0])
    simulated_doublet_score = {}
    if kwargs["batch"] != None:
        for batch in np.unique(config.adata.obs[kwargs["batch"]]):
            sub = config.adata.obs[kwargs["batch"]].values == batch
            subX, sub_doublet_scores_obs_, sub_doublet_scores_sim_ = scrublet_(X[sub,:], kwargs)
            doublet_scores_obs[sub] = sub_doublet_scores_obs_
            X_umap[sub,:] = subX
            simulated_doublet_score[batch] = sub_doublet_scores_sim_
    else:
        subX, sub_doublet_scores_obs_, sub_doublet_scores_sim_ = scrublet_(X, kwargs)
        doublet_scores_obs = sub_doublet_scores_obs_
        X_umap = subX
        simulated_doublet_score[" "] = sub_doublet_scores_sim_

    d = {
        "obsm" : X_umap,
        "obs" : {
            "doublet_scores_obs": doublet_scores_obs,
            "doublet_scores_obs--keep": np.ones(X.shape[0], bool)
        },
        "uns" : {"simulated_scores": simulated_doublet_score}
    }

    return d

def scrublet_reset_lims():

    batch = config.adata.uns[config.selected]["parameters"]["batch"]
    thresholds = config.adata.uns[config.selected]["parameters"]["thresholds"]
    for l in config.adata.uns[config.selected]["parameters"]["thresholds"]:
        name = get_name("doublet_scores_obs")
        name_keep = get_name("doublet_scores_obs")+"--keep"
        if batch != None:
            sub = config.adata.obs[batch] == l["batch"]
            config.adata.obs[name_keep].values[sub]  = config.adata.obs[name].values[sub] <= float(l["max"])
        else:
            config.adata.obs[name_keep]  = config.adata.obs[name].values <= float(l["max"])

def scrublet_plot():

    scrublet_reset_lims()

    pos = get_node_pos(config.selected)
    node = get_node(config.selected)
    if not node["data"]["computed"]:
        return []

    l = []

    for b,c_sim in config.adata.uns[config.selected]["simulated_scores"].items():

        if b == " ":
            sub = np.ones(config.adata.X.shape[0],dtype=bool)
            lims_max = float([i for i in node["data"]["parameters"]["thresholds"]][0]["max"])
            x = config.adata.uns[config.selected]["simulated_scores"][" "]
        else:
            sub = config.adata.obs[node["data"]["parameters"]["batch"]] == b
            lims_max = float([i for i in node["data"]["parameters"]["thresholds"] if i["batch"] == b][0]["max"])
            x = config.adata.uns[config.selected]["simulated_scores"][b]

        # res = int(get_table_value(data,b,"scrublet nBins"))
        X = config.adata.obsm[config.selected][sub,:]
        c = np.array(config.adata.obs[f"{config.selected}--doublet_scores_obs"])[sub]
        if node["data"]["plot"]["show_scores"]:
            c = c
        else:
            c = c > lims_max
        order = np.argsort(c)

        fig = px.histogram(x=np.array(config.adata.obs[f"{config.selected}--doublet_scores_obs"])[sub],histnorm="probability",barmode="overlay")
        fig.add_traces(list(
            px.histogram(x=x,
                         color_discrete_sequence=["red"],
                         histnorm="probability",
                         barmode="overlay"
            ).select_traces()
        ))
        fig.add_traces(list(
            px.line(x=[lims_max, lims_max],
                       y=[0, .1],
                       color_discrete_sequence=["black"],
            ).select_traces()
        ))

        l += [dbc.Row([
                dbc.Col([
                    dcc.Graph(
                        figure= fig,
                        style={"width": "80vh", "height": "60vh"}
                    ),
                ],
                ),
                dbc.Col([
                    dcc.Graph(
                        figure= px.scatter(
                                        x=X[order,0],
                                        y=X[order,1],
                                        color=c[order],
                                    ),
                        style={"width": "60vh", "height": "60vh"}
                    ),
                ],
                ),
            ])
        ]

    return l

config.methods["scrublet"] = {

    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": args,

    "function":scrublet_f,

    "plot":scrublet_plot,

}