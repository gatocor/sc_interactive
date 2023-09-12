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

from app import app

from .. import config

def args_scrublet():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_batch = get_batch_keys()

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
            "description":"str, optional (default: None) Batch key to use. The Doublet metric will be computed independently in each set of cells separated by batch. If None, use the full dataset.",
            "value":None,
            "clearable":True,
            "options":options_batch,
            "summary":True
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
    ]

@app.callback(
    dash.Output("analysis_distance_metric","options", allow_duplicate = True),
    dash.Output("analysis_distance_metric","value", allow_duplicate = True),
    dash.Input("analysis_use_approx_neighbors","on"),
    prevent_initial_call=True
)
def metric_options(approximate):

    if approximate:
        return ["euclidean","manhattan","angular","hamming","dot"], "euclidean"
    else:
        return ["euclidean","cityblock","cosine","haversine","l1","l2","manhattan"], "euclidean"

def scrublet_(X, kwargs):

    scrub = scrublet.Scrublet(X)

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

    return X, scrub.doublet_scores_obs_, scrub.doublet_scores_sim_

def f_scrublet(name_analysis, kwargs):

    pos = get_node_pos(config.selected)

    if kwargs['batch_key'] == None:

        X, doublet_scores_obs_, doublet_scores_sim_ = scrublet_(config.adata.X, kwargs)

        config.adata.obs[name_analysis+"_scrublet_score"] = doublet_scores_obs_
        config.adata.obsm[name_analysis+"_UMAP"] = X
        config.adata.uns["sc_interactive"][name_analysis] = {
            "doublets_simulated_scrublet_score" : doublet_scores_sim_
        }

        config.graph[pos]['data']['batch'] = None

    else:

        config.adata.obs[name_analysis+"_scrublet_score"] = -1
        config.adata.obsm[name_analysis+"_UMAP"] = np.zeros([config.adata.shape[0],2])
        config.adata.uns["sc_interactive"][name_analysis] = {
            "doublets_simulated_scrublet_score" : {}
        }

        for b in config.adata.obs[kwargs['batch_key']].unique():

            sub = config.adata.obs[kwargs['batch_key']].values == b

            X, doublet_scores_obs_, doublet_scores_sim_ = scrublet_(config.adata.X[sub,:], kwargs)

            config.adata.obs.loc[sub,name_analysis+"_scrublet_score"] = doublet_scores_obs_
            config.adata.obsm[name_analysis+"_UMAP"][sub,:] = X
            config.adata.uns["sc_interactive"][name_analysis]["doublets_simulated_scrublet_score"][b] = doublet_scores_sim_

        config.graph[pos]['data']['batch'] = config.adata.obs[kwargs['batch_key']].values

    #Make empty table
    add = ["max","nBins"]
    if kwargs["batch_key"] == None:
        rows = [None]
    else:
        rows = np.sort(config.adata.obs[kwargs["batch_key"]].unique())
    
    columns, data = make_thresholds_table(['scrublet'], rows, add)

    #Fill table 
    for i in rows:
        set_table_value(data, i, "scrublet max", 1)
        set_table_value(data, i, "scrublet nBins", 20)

    config.graph[pos]['data']['doublets_simulated_scrublet_score'] = config.adata.uns["sc_interactive"][name_analysis]["doublets_simulated_scrublet_score"].copy()
    config.graph[pos]['data']['doublets_score'] = config.adata.obs[name_analysis+"_scrublet_score"].values
    config.graph[pos]['data']['UMAP'] = config.adata.obsm[name_analysis+"_UMAP"].copy()
    config.graph[pos]['data']['threshold'] = {'columns':columns,'data':data}
    config.graph[pos]['data']['filter'] = np.ones_like(config.adata.X.shape[0])>0
    config.graph[pos]['data']['show_scores'] = True

def rm_scrublet(name_analysis):

    config.adata.obs.drop(name_analysis+"_scrublet_score", axis=1, inplace=True)
    del config.adata.obsm[name_analysis+"_UMAP"]
    del config.adata.uns[name_analysis]

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

    node = get_node(name_analysis)
    if not node['data']['computed']:
        return []

    node = get_node(name_analysis)
    pos = get_node_pos(name_analysis)

    columns = node['data']['threshold']['columns']
    data = node['data']['threshold']['data']

    l = [
        html.H1("Scrublet analysis"),
        dbc.Row(
            dash_table.DataTable(
                id="scrublet_threshold_table",
                columns=columns,
                data=data
            )
        ),
        dbc.Row(
            [
                dbc.Col(),
                dbc.Col(
                    daq.BooleanSwitch(id="scrublet_toggle",label="Removed Cells/Scrublet Score",on=node['data']['show_scores']),
                )
            ]
        )
    ]

    if node['data']['parameters']['batch_key'] == None:
        lims_max = float(get_table_column(data,"scrublet max")[0])
        res = int(get_table_column(data,"scrublet nBins")[0])
        X = np.array(config.graph[pos]['data']['UMAP'])
        c = np.array(config.graph[pos]['data']['doublets_score'])
        if node['data']['show_scores']:
            c = c
        else:
            c = c > lims_max
        order = np.argsort(c)
        c_sim = config.graph[pos]['data']['doublets_simulated_scrublet_score']
        plot_max = hist_vline(c_sim, res)
        l += [
                dbc.Row([
                    dbc.Col([
                        dcc.Graph(
                            figure={
                                    "data":[
                                        go.Histogram(
                                            x=c_sim,
                                            nbinsx=res,
                                            name='simulated doublets',
                                        ),
                                        go.Scatter(
                                            x=[lims_max, lims_max],
                                            y=[0,plot_max],
                                            name='Min threshold',
                                            marker=dict(color='orange'),
                                            opacity=0.7
                                        ),
                                    ],
                                    "layout":{
                                            'xlabel':'scrublet_score'
                                            # 'yaxis':{
                                            #     'scaleanchor':"x",
                                            #     'scaleratio':1,
                                            # },
                                    }
                                },
                            style={'width': '90vh', 'height': '60vh'}
                        ),
                    ],
                    ),
                    dbc.Col([
                        dcc.Graph(
                            id="scrublet_scatter",
                            figure={
                                    "data":[
                                        go.Scatter(
                                            x=X[order,0],
                                            y=X[order,1],
                                            marker={'color':qualitative_colors(c[order])},
                                            mode='markers',
                                        )
                                    ],
                                    "layout":{
                                            'yaxis':{
                                                'scaleanchor':"x",
                                                'scaleratio':1,
                                            },
                                    }
                                },
                            style={'width': '60vh', 'height': '60vh'}
                        ),
                    ],
                    ),
                ],
                justify='center'
                )
            ]
        
    else:
        for b,c_sim in config.graph[pos]['data']['doublets_simulated_scrublet_score'].items():
            sub = np.array(config.graph[pos]['data']['batch']) == b
            lims_max = float(get_table_value(data,b,"scrublet max"))
            res = int(get_table_value(data,b,"scrublet nBins"))
            X = np.array(config.graph[pos]['data']["UMAP"])[sub,:]
            c = np.array(config.graph[pos]['data']["doublets_score"])[sub]
            if node['data']['show_scores']:
                c = c
            else:
                c = c > lims_max
            order = np.argsort(c)

            plot_max = hist_vline(c_sim, res)
            l += [
                    dbc.Row([
                        dbc.Col([
                            dcc.Graph(
                                figure={
                                        "data":[
                                            go.Histogram(
                                                x=c_sim,
                                                name='Simulated Doublets',
                                                nbinsx=res
                                            ),
                                            go.Scatter(
                                                x=[lims_max, lims_max],
                                                y=[0,plot_max],
                                                name='Max threshold',
                                                marker=dict(color='red'),
                                                opacity=0.7
                                            ),
                                        ],
                                        "layout":{
                                                'xlabel':'scrublet_score'
                                                # 'yaxis':{
                                                #     'scaleanchor':"x",
                                                #     'scaleratio':1,
                                                # },
                                        }
                                    },
                                style={'width': '90vh', 'height': '60vh'}
                            ),
                        ],
                        align='center'
                        ),
                        dbc.Col([
                            dcc.Graph(
                                figure={
                                        "data":[
                                            go.Scatter(
                                                x=X[order,0],
                                                y=X[order,1],
                                                marker={'color':qualitative_colors(c[order])},
                                                mode='markers',
                                            )
                                        ],
                                        "layout":{
                                                'yaxis':{
                                                    'scaleanchor':"x",
                                                    'scaleratio':1,
                                                },
                                        }
                                    },
                                style={'width': '60vh', 'height': '60vh'}
                            ),
                        ],
                        ),
                    ],
                    justify='center'
                    )
                ]

    return l

@app.callback(
    dash.Output("analysis_plot","children",allow_duplicate=True),
    dash.Input("scrublet_threshold_table","data"),
    prevent_initial_call=True
)
def update_table(data):

    prevent_race('scrublet')

    pos = get_node_pos(config.selected)
    config.graph[pos]['data']['threshold']['data'] = data
    batch = get_node_parameters(config.selected)['batch_key']

    if batch == None:

        col = np.array([float(i) for i in get_table_column(data,"scrublet max")])
        s = np.array(config.graph[pos]['data']['doublets_score']) > float(get_table_value(data,batch,"scrublet max"))
        s <= 1

        config.graph[pos]['data']['filter'] = s

    elif batch != None:

        col = np.array([float(i) for i in get_table_column(data,"scrublet max")])
        s = np.ones_like(config.graph[pos]['data']['doublets_score'])
        for b in np.unique(config.graph[pos]['data']['batch']):    
            sub = np.array(config.graph[pos]['data']['batch']) == b
            s[sub] = np.array(config.graph[pos]['data']['doublets_score'])[sub] <= float(get_table_value(data,b,"scrublet max"))
            s <= 1

        config.graph[pos]['data']['filter'] = s

    return plot_scrublet(config.selected)


@app.callback(
    dash.Output("analysis_plot","children",allow_duplicate=True),
    dash.Output("scrublet_toggle","on",allow_duplicate=True),
    dash.Input("scrublet_toggle","on"),
    prevent_initial_call=True
)
def update_table(data):
    
    prevent_race('scrublet')

    pos = get_node_pos(config.selected)
    config.graph[pos]['data']['show_scores'] = data

    return plot_scrublet(config.selected), data
