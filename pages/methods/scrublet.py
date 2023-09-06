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

from ..functions import *

from app import app

from .. import config

def args_scrublet():

    options = node_names(exclude_downstream_from_node=config.selected) 
    
    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options":options
        },
        {
            "input":"Dropdown",
            "name":"batch_key",
            "description":"str, optional (default: 'Full') Batch key to use. The Doublet metric will be computed independently in each set of cells separated by batch. If None, use the full dataset.",
            "value":'Full',
            "clearable":False,
            "options":options
        },
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
            "value":"correlation",
            "clearable":False,
            "options":["euclidean","manhattan","correlation"]
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
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"svd_solver",
            "description":"str, optional (default: 'arpack') SVD solver to use. See available options for  'svd_solver' from 'sklearn.decomposition.PCA' or 'algorithm' from 'sklearn.decomposition.TruncatedSVD'",
            "value": "arpack",
            "clearable":False,
            "options": ["arpack"]
        }
    ]

def f_scrublet(name_analysis, kwargs):

    scrub = scrublet.Scrublet(config.adata.X)

    scrub.scrub_doublets(
        synthetic_doublet_umi_subsampling = kwargs["synthetic_doublet_umi_subsampling"], 
        use_approx_neighbors = kwargs["use_approx_neighbors"], 
        distance_metric = kwargs["distance_metric"], 
        min_counts = kwargs["min_counts"],
        min_cells = kwargs["min_cells"], 
        min_gene_variability_pctl = kwargs["min_gene_variability_pctl"], 
        log_transform = kwargs["log_transform"], 
        mean_center = kwargs["mean_center"], 
        normalize_variance = kwargs["normalize_variance"], 
        n_prin_comps = kwargs["n_prin_comps"], 
        svd_solver= kwargs["svd_solver"]
    )
            
    X = scrublet.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)

    config.adata.obs[name_analysis+"_scrublet_score"] = scrub.doublet_scores_obs_
    config.adata.obsm[name_analysis+"_UMAP"] = X
    config.adata.uns["sc_interactive"][name_analysis] = {
        "doublets_simulated_scrublet_score" : scrub.doublet_scores_sim_
    }

def rm_scrublet(name_analysis):

    config.adata.obs.drop(name_analysis+"_scrublet_score", axis=1, inplace=True)
    del config.adata.obsm[name_analysis+"_UMAP"]
    del config.adata.uns["sc_interactive"][name_analysis]

def rename_scrublet(name_analysis, name_new_analysis):

    ls = [name_analysis+"_scrublet_score"]
    cols = {}    
    for i in config.adata.obs.columns.values:
        if i not in ls:
            cols[i] = i
        else:
            name = i.split(name_analysis)[-1]
            name = (name_new_analysis+"_"+name)
            cols[i] = name

    config.adata.obs.rename(columns=cols, inplace=True)

def plot_scrublet(name_analysis):

    l = []
    # dic = get_node_parameters(name_analysis,str2list=True)['measure'] 

    # for metric in dic:

    #     var_selected_data = name_analysis+"_"+metric
            
    #     l += [
    #             dbc.Row(
    #                 [dcc.Graph(id="Histogram",
    #                         figure={
    #                                 "data":[
    #                                     go.Histogram(
    #                                         x=config.adata.obs[var_selected_data].values,
    #                                         nbinsx=100,
    #                                         name='Histogram',
    #                                         marker=dict(color='blue'),
    #                                         opacity=0.7
    #                                     ),
    #                                 ],
    #                                 "layout":{
    #                                         'title': f'Histogram of {var_selected_data}',
    #                                         'xaxis': {'title': var_selected_data},
    #                                         'yaxis': {'title': 'Count'},
    #                                         'barmode': 'overlay',
    #                                         'width':1500,
    #                                         'height':400,
    #                                 }
    #                             }
    #                 )],
    #                 justify="center",
    #                 style={'width': '90%', 'margin': 'auto'}
    #             )
    #         ]
        
    return l
