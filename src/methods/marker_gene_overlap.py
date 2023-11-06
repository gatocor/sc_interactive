import scanpy as sc

from general import *

marker_gene_overlap_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"Dropdown",
            "name":"var_key",
            "description":"A columns of adata.var",
            "properties" : {
                "value":{"function":"[i for i in config.adata.var.columns.values][0]"},
                "options":{"function":"[i for i in config.adata.var.columns.values]"},
            }
        },
        {
            "input":"AgTable",
            "name":"reference_markers",
            "description":"A marker gene dictionary object. Keys should be strings with the cell identity name and values are sets or lists of strings which match format of var_key.",
            "properties":{
                "header":[
                    { "headerName": "Type", "field":"type", "editable": True },
                    { "headerName": "Markers", "field":"markers", "editable": True },
                ],
                "data":[{"type":"example_cell_type","markers":"T,Sox2"}],
            },
            "addRows":{"type":"", "field":""},
            "deleteRows": True                
        },
        {
            "input":"Dropdown",
            "name":"method",
            "description":"Method to calculate marker gene overlap. 'overlap_count' uses the intersection of the gene set, 'overlap_coef' uses the overlap coefficient, and 'jaccard' uses the Jaccard index.",
            "properties" : {
                "value":"overlap_count",
                "clearable":False,
                "options":['overlap_count', 'overlap_coef', 'jaccard'],
            }
        },
        {
            "input":"Dropdown",
            "name":"normalize",
            "description":"Normalization option for the marker gene overlap output. This parameter can only be set when method is set to 'overlap_count'. 'reference' normalizes the data by the total number of marker genes given in the reference annotation per group. 'data' normalizes the data by the total number of marker genes used for each cluster.",
            "properties" : {
                "value":None,
                "clearable":True,
                "options":['reference', 'data'],
            }
        },
        {
            "input":"Input",
            "name":"top_n_markers",
            "description":"The number of top data-derived marker genes to use. By default the top 100 marker genes are used. If adj_pval_threshold is set along with top_n_markers, then adj_pval_threshold is ignored.",
            "properties" : {
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"adj_pval_threshold",
            "description":"A significance threshold on the adjusted p-values to select marker genes. This can only be used when adjusted p-values are calculated by sc.tl.rank_genes_groups(). If adj_pval_threshold is set along with top_n_markers, then adj_pval_threshold is ignored.",
            "properties" : {
                "value":None,
                "type":"number"
            }
        },
    ],

    "postexecution" : [],

    "plot" : [
        # {
        #     "input":"Dropdown",
        #     "name":"plot_style",
        #     "description":"Style of plot.",
        #     "properties":{
        #         "value":"heatmap",
        #         "clearable":False,
        #         "options":["heatmap","scattermap","table"]
        #     },
        # },
        # {
        #     "input":"Dropdown",
        #     "name":"plot_n_genes",
        #     "description":"Number of genes ploted.",
        #     "properties":{
        #         "value":2,
        #         "clearable":False,
        #         "options":np.arange(1,100)
        #     },
        # },
        # {
        #     "input":"Dropdown",
        #     "name":"var_key",
        #     "description":"var key to show.",
        #     "properties":{
        #         "value":{"function":"config.adata.var.columns.values[0]"},
        #         "clearable":False,
        #         "options":{"function":"config.adata.var.columns.values"}
        #     },
        # },
        # {
        #     "input":"Dropdown",
        #     "name":"values_to_plot",
        #     "description":"Values to plot.",
        #     "properties":{
        #         "value":"scores",
        #         "clearable":False,
        #         "options":['scores', 'logfoldchanges', 'pvals', 'pvals_adj']
        #     },
        # },
        # {
        #     "input":"Dropdown",
        #     "name":"plot_cluster",
        #     "description":"Values to plot.",
        #     "properties":{
        #         "value":{"function":"config.adata.uns['marker_gene_overlap']['scores'].dtype.names[0]"},
        #         "clearable":False,
        #         "options":{"function":"config.adata.uns['marker_gene_overlap']['scores'].dtype.names"}
        #     },
        # },
    ]
}

def marker_gene_overlap_f(adata, kwargs):

    d = pd.DataFrame()
    d["source"] = config.adata.var[kwargs["var_key"]].values
    d["target"] = config.adata.var.index.values
    d.set_index("source",inplace=True)
    
    dic = {}
    for i in kwargs["reference_markers"]:
        dic[i["type"]] = {d.loc[j,"target"] for j in i["markers"].split(",")}

    d = sc.tl.marker_gene_overlap(
                config.adata,
                reference_markers=dic,
                method=kwargs["method"],
                normalize=kwargs["normalize"],
                top_n_markers=kwargs["top_n_markers"],
                adj_pval_threshold=kwargs["adj_pval_threshold"],
                inplace=False,
                )
    
    config.adata.uns["marker_gene_overlap"] = d

def marker_gene_overlap_plot():

    fig = px.imshow(config.adata.uns["marker_gene_overlap"].values, 
                    y=config.adata.uns["marker_gene_overlap"].index.values,
                    x=config.adata.uns["marker_gene_overlap"].columns.values,
                    text_auto=True
                    )

    fig = dcc.Graph(figure=fig)

    return fig

config.methods["marker_gene_overlap"] = {

    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": marker_gene_overlap_args,

    "function":marker_gene_overlap_f,

    "plot":marker_gene_overlap_plot,

}