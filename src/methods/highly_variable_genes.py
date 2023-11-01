import scanpy as sc
import dash_bootstrap_components as dbc
import plotly.express as px

from general import *

highly_variable_genes_args = {
    
    "execution" : [
        ARGINPUT,
        {
            "input":"Input",
            "name":"n_top_genes",
            "description":"Number of highly-variable genes to keep. Mandatory if flavor='seurat_v3'.",
            "properties":{
                "value":1000,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"min_disp",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "properties":{
                "value":.5,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"max_disp",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "properties":{
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"min_mean",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "properties":{
                "value":0.0125,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"max_mean",
            "description":"If n_top_genes unequals None, this and all other cutoffs for the means and the normalized dispersions are ignored. Ignored if flavor='seurat_v3'.",
            "properties":{
                "value":3,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"span",
            "description":"The fraction of the data (cells) used when estimating the variance in the loess model fit if flavor='seurat_v3'.",
            "properties":{
                "value":.3,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"n_bins",
            "description":"Number of bins for binning the mean gene expression. Normalization is done with respect to each bin. If just a single gene falls into a bin, the normalized dispersion is artificially set to 1. You’ll be informed about this if you set settings.verbosity = 4.",
            "properties":{
                "value":20,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"flavor",
            "description":"Choose the flavor for identifying highly variable genes. For the dispersion based methods in their default workflows, Seurat passes the cutoffs whereas Cell Ranger passes n_top_genes.",
            "properties":{
                "value":"seurat",
                "clearable":False,
                "options":["seurat", "cell_ranger", "seurat_v3"] 
            }
        },
        {
            "input":"Dropdown",
            "name":"batch_key",
            "description":"If specified, highly-variable genes are selected within each batch separately and merged. This simple process avoids the selection of batch-specific genes and acts as a lightweight batch correction method. For all flavors, genes are first sorted by how many batches they are a HVG. For dispersion-based flavors ties are broken by normalized dispersion. If flavor = 'seurat_v3', ties are broken by the median (across batches) rank based on within-batch normalized variance.",
            "properties":{
                "value":None,
                "clearable":True,
                "options":{"function":"[i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in ['int','bool','category','object','str']]"} 
            }
        },
    ],

    "postexecution" : [],

    "plot" : []
}

def highly_variable_genes_f(adata, kwargs):

    sc.pp.highly_variable_genes(config.adata,
        n_top_genes = kwargs["n_top_genes"], 
        min_disp = kwargs["min_disp"], 
        max_disp = kwargs["max_disp"], 
        min_mean = kwargs["min_mean"], 
        max_mean = kwargs["max_mean"], 
        span = kwargs["span"], 
        n_bins = kwargs["n_bins"], 
        flavor = kwargs["flavor"], 
    )

def highly_variable_genes_plot():

    m = config.adata.var["means"]
    c = config.adata.var["highly_variable"]
    if "dispersions" in config.adata.var.keys():
        v = config.adata.var["dispersions"]
        vn = config.adata.var["dispersions_norm"]
        p = "Dispersions"
    else:
        v = config.adata.var["variances"]
        vn = config.adata.var["variances_norm"]
        p = "Variances"

    fig = px.scatter(x=m, y=v, color=c)    
    fig.layout["xaxis"]["title"] = "Means"
    fig.layout["yaxis"]["title"] = p
    fig.layout["legend"]["title"] = "Selected genes"

    fign = px.scatter(x=m, y=vn, color=c)    
    fign.layout["xaxis"]["title"] = "Means"
    fign.layout["yaxis"]["title"] = p
    fign.layout["legend"]["title"] = "Selected genes"

    l = dbc.Row([
            dbc.Col(
                plot_center(
                    dcc.Graph(
                        figure = fig
                    )
                ),
            ),
            dbc.Col(
                plot_center(
                    dcc.Graph(
                        figure = fign
                    )
                ),
            ),
        ])   
                
    return l

config.methods["highly_variable_genes"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": highly_variable_genes_args,

    "function": highly_variable_genes_f,

    "plot": highly_variable_genes_plot,

}