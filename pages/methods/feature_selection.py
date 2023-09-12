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

def args_feature_selection():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_batch = get_batch_keys()

    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Representation to use as input of the method.",
            "value":None,
            "clearable":True,
            "options":options
        },
        {
            "input":"Dropdown",
            "name":"batch_key",
            "description":"If specified, highly-variable genes are selected within each batch separately and merged. This simple process avoids the selection of batch-specific genes and acts as a lightweight batch correction method. For all flavors, genes are first sorted by how many batches they are a HVG. For dispersion-based flavors ties are broken by normalized dispersion. If flavor = 'seurat_v3', ties are broken by the median (across batches) rank based on within-batch normalized variance.",
            "value":None,
            "clearable":True,
            "options":options_batch
        },
        {
            "input":"Input",
            "name":"n_top_genes",
            "description":"Number of highly-variable genes to keep. Mandatory if flavor='seurat_v3'.",
            "value":1000,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"min_disp",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "value":.5,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"max_disp",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "value":None,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"min_mean",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "value":0.0125,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"max_mean",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "value":3,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"span",
            "description":"The fraction of the data (cells) used when estimating the variance in the loess model fit if flavor='seurat_v3'.",
            "value":.3,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"n_bins",
            "description":"Number of bins for binning the mean gene expression. Normalization is done with respect to each bin. If just a single gene falls into a bin, the normalized dispersion is artificially set to 1. You’ll be informed about this if you set settings.verbosity = 4.",
            "value":20,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"flavor",
            "description":"Choose the flavor for identifying highly variable genes. For the dispersion based methods in their default workflows, Seurat passes the cutoffs whereas Cell Ranger passes n_top_genes.",
            "value":'seurat',
            "clearable":False,
            "options":['seurat', 'cell_ranger', 'seurat_v3'] 
        },
    ]

def f_feature_selection(name_analysis, kwargs):
        
    adata_copy = sc.AnnData(X=config.adata.X.copy())

    if kwargs["flavor"] != "seurat_v3":
        sc.pp.log1p(adata_copy)

    kwargs_copy = kwargs.copy()
    del kwargs_copy["input"]
    sc.pp.highly_variable_genes(adata_copy,
        **kwargs_copy
    )
    
    config.adata.var[name_analysis+"_highly_variable"] = adata_copy.var["highly_variable"].values
    config.adata.var[name_analysis+"_means"] = adata_copy.var["means"].values

    if "dispersions" in adata_copy.var.columns:
        config.adata.var[name_analysis+"_dispersions"] = adata_copy.var["dispersions"].values
        config.adata.var[name_analysis+"_dispersions_norm"] = adata_copy.var["dispersions_norm"].values
    else:
        config.adata.var[name_analysis+"_variances"] = adata_copy.var["variances"].values
        config.adata.var[name_analysis+"_variances_norm"] = adata_copy.var["variances_norm"].values

def rm_feature_selection(name_analysis):

    config.adata.var.drop(name_analysis+"_highly_variable", axis=1)
    config.adata.var.drop(name_analysis+"_means", axis=1)

    if name_analysis+"_dispersions" in config.adata.var.columns:
        config.adata.var.drop(name_analysis+"_dispersions", axis=1)
        config.adata.var.drop(name_analysis+"_dispersions_norm", axis=1)
    else:
        config.adata.var.drop(name_analysis+"_variances",axis=1)
        config.adata.var.drop(name_analysis+"_variances_norm",axis=1)

    return

def rename_feature_selection(name_analysis, name_new_analysis):

    config.adata.var[name_new_analysis+"_highly_variable"] = config.adata.var[name_analysis+"_highly_variable"]
    config.adata.var[name_new_analysis+"_means"] = config.adata.var[name_analysis+"_means"]
    config.adata.var.drop(name_analysis+"_highly_variable", axis=1)
    config.adata.var.drop(name_analysis+"_means", axis=1)

    if "dispersions" in config.adata.var.columns:
        config.adata.var[name_new_analysis+"_dispersions"] = config.adata.var[name_analysis+"_dispersions"]
        config.adata.var[name_new_analysis+"_dispersions_norm"] = config.adata.var[name_analysis+"_dispersions_norm"]
        config.adata.var.drop(name_analysis+"_dispersions", axis=1)
        config.adata.var.drop(name_analysis+"_dispersions_norm", axis=1)
    else:
        config.adata.var[name_new_analysis+"_variances"] = config.adata.var[name_analysis+"_variances"]
        config.adata.var[name_new_analysis+"_variances_norm"] = config.adata.var[name_analysis+"_variances_norm"]
        config.adata.var.drop(name_analysis+"_variances",axis=1)
        config.adata.var.drop(name_analysis+"_variances_norm",axis=1)

    return

def plot_feature_selection(name_analysis):

    node = get_node(config.selected)
    if node['data']['computed']:

        m = config.adata.var[name_analysis+"_means"].values
        c =config.adata.var[name_analysis+"_highly_variable"].values
        if name_analysis+"_dispersions" in config.adata.var.columns:
            v = config.adata.var[name_analysis+"_dispersions"].values
            vn = config.adata.var[name_analysis+"_dispersions_norm"].values
        else:
            v = config.adata.var[name_analysis+"_variances"].values
            vn = config.adata.var[name_analysis+"_variances_norm"].values

        color_map = {
            True:"orange",
            False:"blue"
        }

        return [
            dbc.Col(
                dcc.Graph(
                    figure = {'data':[
                                go.Scattergl(
                                    x=m,
                                    y=v,
                                    mode='markers',
                                    name='Min threshold',
                                    marker=dict(
                                        color=[color_map[i] for i in c],
                                    ),
                                )],
                            #  'layout':{
                            #          'title': "Unnormalized analysis",
                            #          'xaxis': "Mean",
                            #          'yaxis': "Dispersion",
                            #          'barmode': 'overlay'
                            #  }
                            },
                )
            ),
            dbc.Col(
                dcc.Graph(
                    figure = {'data':[
                                go.Scattergl(
                                    x=m,
                                    y=vn,
                                    mode='markers',
                                    name='Min threshold',
                                    marker=dict(
                                        color=[color_map[i] for i in c],
                                    ),
                                )],
                            #   'layout':{
                            #           'title': "Normalized analysis",
                            #           'xaxis': "Mean",
                            #           'yaxis': "Normalized Dispersion",
                            #   }
                            },
                )
            ),
        ]
    
    else:

        return []