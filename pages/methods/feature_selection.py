import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import plotly.graph_objs as go

def feature_selection_args(adata):

    options = []
    if "__interactive__" in adata.uns.keys():
        options = ["Raw"]+[i[2:] for i in adata.obsm.keys()] 

    return [
        "Feature Selection",
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Representation to use as input of the method.",
            "value":"Raw",
            "clearable":True,
            "options":options
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
            "value":100000000000000,
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
        {
            "input":"Dropdown",
            "name":"batch_key",
            "description":"If specified, highly-variable genes are selected within each batch separately and merged. This simple process avoids the selection of batch-specific genes and acts as a lightweight batch correction method. For all flavors, genes are first sorted by how many batches they are a HVG. For dispersion-based flavors ties are broken by normalized dispersion. If flavor = 'seurat_v3', ties are broken by the median (across batches) rank based on within-batch normalized variance.",
            "value":None,
            "clearable":True,
            "options":[str(i) for i in adata.obs.columns.values if (adata.obs.dtypes[i] in ["category" ,object, str, int])]
        },
    ]

def f_feature_selection(adata, name_analysis, **kwargs):
        
    if kwargs["input"] == "Raw":
        adata_copy = sc.AnnData(X=adata.X.copy())
    else:
        adata_copy = sc.AnnData(X=adata.obsm["X_"+kwargs["input"]].copy())

    if kwargs["flavor"] != "seurat_v3":
        sc.pp.log1p(adata_copy)

    kwargs_copy = kwargs.copy()
    del kwargs_copy["input"]
    sc.pp.highly_variable_genes(adata_copy,
        **kwargs_copy
    )
    
    adata.var[name_analysis+"_highly_variable"] = adata_copy.var["highly_variable"].values
    adata.var[name_analysis+"_means"] = adata_copy.var["means"].values

    if "dispersions" in adata_copy.var.columns:
        adata.var[name_analysis+"_dispersions"] = adata_copy.var["dispersions"].values
        adata.var[name_analysis+"_dispersions_norm"] = adata_copy.var["dispersions_norm"].values
    else:
        adata.var[name_analysis+"_variances"] = adata_copy.var["variances"].values
        adata.var[name_analysis+"_variances_norm"] = adata_copy.var["variances_norm"].values

    adata.obsm["X_"+name_analysis] = adata_copy.X[:,adata_copy.var["highly_variable"]]

    adata.uns[name_analysis] = adata_copy.uns["hvg"]

def make_feature_selection_plots1(adata, name_analysis):

    if name_analysis+"_means" in adata.var.columns.values:

        m = adata.var[name_analysis+"_means"].values
        c =adata.var[name_analysis+"_highly_variable"].values
        if name_analysis+"_dispersions" in adata.var.columns:
            v = adata.var[name_analysis+"_dispersions"].values
            vn = adata.var[name_analysis+"_dispersions_norm"].values
        else:
            v = adata.var[name_analysis+"_variances"].values
            vn = adata.var[name_analysis+"_variances_norm"].values

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
                                )]
                            }
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
                                )]
                            }
                )
            ),
        ]
    else:
        return []

def make_feature_selection_plots2(adata, name_analysis):
    return []
