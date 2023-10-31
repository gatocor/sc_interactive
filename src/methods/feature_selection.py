import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
from dash import dcc, dash_table
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash

from general import *

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
            "description":"str, optional (default: None) Batch key to use. The highly varying will be computed independently in each set of cells separated by batch. If None, use the full dataset.",
            "value":None,
            "clearable":True,
            "options":options_batch,
            "summary":True
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
            "value":"seurat",
            "clearable":False,
            "options":["seurat", "cell_ranger", "seurat_v3"] 
        },
    ]

def f_feature_selection(name_analysis, kwargs, sub_name, sub):

    l = [i["data"]["method"] for i in get_node_ancestors(name_analysis)]

    adata_copy = sc.AnnData(X=config.adata.X[sub,:])

    if kwargs["flavor"] == "seurat_v3" and "log1p" in l:
        np.expm1(adata_copy)
        del adata_copy.uns["log1p"]
    elif "log1p" not in l:
        sc.pp.log1p(adata_copy)

    kwargs_copy = kwargs.copy()
    del kwargs_copy["input"]
    del kwargs_copy["batch_key"]
    sc.pp.highly_variable_genes(adata_copy,
        **kwargs_copy,
    )

    pos = get_node_pos(name_analysis)
    if "var" not in config.graph[pos]["data"].keys():
        config.graph[pos]["data"]["var"] = {}

    config.graph[pos]["data"]["var"]["means_"+str(sub_name)] = adata_copy.var["means"].values
    config.graph[pos]["data"]["var"]["highly_variable_"+str(sub_name)] = adata_copy.var["highly_variable"].values
    if "dispersions" in adata_copy.var.columns:
        config.graph[pos]["data"]["var"]["dispersions_"+str(sub_name)] = adata_copy.var["dispersions"].values
        config.graph[pos]["data"]["var"]["dispersions_norm_"+str(sub_name)] = adata_copy.var["dispersions_norm"].values
    else:
        config.graph[pos]["data"]["var"]["variances_"+str(sub_name)] = adata_copy.var["variances"].values
        config.graph[pos]["data"]["var"]["variances_norm_"+str(sub_name)] = adata_copy.var["variances_norm"].values

def rm_feature_selection(name_analysis):

    # batch = get_node(name_analysis)["data"]["parameters"]["batch_key"]

    # if batch == None:
    #     config.adata.var.drop(name_analysis+"_highly_variable", axis=1)
    #     config.adata.var.drop(name_analysis+"_means", axis=1)

    #     if name_analysis+"_dispersions" in config.adata.var.columns:
    #         config.adata.var.drop(name_analysis+"_dispersions", axis=1)
    #         config.adata.var.drop(name_analysis+"_dispersions_norm", axis=1)
    #     else:
    #         config.adata.var.drop(name_analysis+"_variances",axis=1)
    #         config.adata.var.drop(name_analysis+"_variances_norm",axis=1)

    return

def rename_feature_selection(name_analysis, name_new_analysis):

    # config.adata.var[name_new_analysis+"_highly_variable"] = config.adata.var[name_analysis+"_highly_variable"]
    # config.adata.var[name_new_analysis+"_means"] = config.adata.var[name_analysis+"_means"]
    # config.adata.var.drop(name_analysis+"_highly_variable", axis=1)
    # config.adata.var.drop(name_analysis+"_means", axis=1)

    # if "dispersions" in config.adata.var.columns:
    #     config.adata.var[name_new_analysis+"_dispersions"] = config.adata.var[name_analysis+"_dispersions"]
    #     config.adata.var[name_new_analysis+"_dispersions_norm"] = config.adata.var[name_analysis+"_dispersions_norm"]
    #     config.adata.var.drop(name_analysis+"_dispersions", axis=1)
    #     config.adata.var.drop(name_analysis+"_dispersions_norm", axis=1)
    # else:
    #     config.adata.var[name_new_analysis+"_variances"] = config.adata.var[name_analysis+"_variances"]
    #     config.adata.var[name_new_analysis+"_variances_norm"] = config.adata.var[name_analysis+"_variances_norm"]
    #     config.adata.var.drop(name_analysis+"_variances",axis=1)
    #     config.adata.var.drop(name_analysis+"_variances_norm",axis=1)

    return

def plot_feature_selection(name_analysis):

    node = get_node(config.selected)
    if not node["data"]["computed"]:
        return []

    l = []

    for i in [i.split("means_")[-1] for i in node["data"]["var"].keys() if "means_" in i]:
        m = node["data"]["var"]["means_"+str(i)]
        c = node["data"]["var"]["highly_variable_"+str(i)]
        if "dispersions_"+str(i) in node["data"]["var"].keys():
            v = node["data"]["var"]["dispersions_"+str(i)]
            vn = node["data"]["var"]["dispersions_norm_"+str(i)]
        else:
            v = node["data"]["var"]["variances_"+str(i)]
            vn = node["data"]["var"]["variances_norm_"+str(i)]

        color_map = {
            True:"orange",
            False:"blue"
        }

        l += [
            dbc.Col(
                dcc.Graph(
                    figure = {"data":[
                                go.Scattergl(
                                    x=m,
                                    y=v,
                                    mode="markers",
                                    name="Min threshold",
                                    marker=dict(
                                        color=[color_map[i] for i in c],
                                    ),
                                )],
                            #  "layout":{
                            #          "title": "Unnormalized analysis",
                            #          "xaxis": "Mean",
                            #          "yaxis": "Dispersion",
                            #          "barmode": "overlay"
                            #  }
                            },
                )
            ),
            dbc.Col(
                dcc.Graph(
                    figure = {"data":[
                                go.Scattergl(
                                    x=m,
                                    y=vn,
                                    mode="markers",
                                    name="Min threshold",
                                    marker=dict(
                                        color=[color_map[i] for i in c],
                                    ),
                                )],
                            #   "layout":{
                            #           "title": "Normalized analysis",
                            #           "xaxis": "Mean",
                            #           "yaxis": "Normalized Dispersion",
                            #   }
                            },
                )
            ),
        ]
            
    return l