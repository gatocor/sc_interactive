import scanpy as sc

from general import *

rank_genes_groups_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"BooleanSwitch",
            "name":"use_raw",
            "description":"Use raw attribute of adata if present.",
            "properties" : {
                "on":True,
            }
        },
        {
            "input":"Dropdown",
            "name":"layer",
            "description":"Key from adata.layers whose value will be used to perform tests on.",
            "properties" : {
                "value":None,
                "clearable":True,
                "options":{"function":"[i for i in config.adata.layers.keys()]"},
            }
        },
        {
            "input":"Input",
            "name":"n_genes",
            "description":"The number of genes that appear in the returned tables. Defaults to all genes.",
            "properties" : {
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"method",
            "description":"The default method is 't-test', 't-test_overestim_var' overestimates variance of each group, 'wilcoxon' uses Wilcoxon rank-sum, 'logreg' uses logistic regression. See [Ntranos18], here and here, for why this is meaningful.",
            "properties" : {
                "value":"t-test",
                "clearable":False,
                "options":["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
            }
        },
        {
            "input":"Dropdown",
            "name":"corr_method",
            "description":"p-value correction method. Used only for 't-test', 't-test_overestim_var', and 'wilcoxon'.",
            "properties" : {
                "value":"benjamini-hochberg",
                "clearable":False,
                "options":["benjamini-hochberg", "bonferroni"],
            }
        },
        {
            "input":"BooleanSwitch",
            "name":"tie_correct",
            "description":"Use tie correction for 'wilcoxon' scores. Used only for 'wilcoxon'.",
            "properties" : {
                "on":False,
            }
        },
        {
            "input":"BooleanSwitch",
            "name":"pts",
            "description":"Compute the fraction of cells expressing the genes.",
            "properties" : {
                "on":False,
            }
        },
    ],

    "postexecution" : [],

    "plot" : [
        {
            "input":"Dropdown",
            "name":"plot_style",
            "description":"Style of plot.",
            "properties":{
                "value":"heatmap",
                "clearable":False,
                "options":["heatmap","scattermap","table"]
            },
        },
        {
            "input":"Dropdown",
            "name":"plot_n_genes",
            "description":"Number of genes ploted.",
            "properties":{
                "value":2,
                "clearable":False,
                "options":np.arange(1,100)
            },
        },
        {
            "input":"Dropdown",
            "name":"var_key",
            "description":"var key to show.",
            "properties":{
                "value":{"function":"config.adata.var.columns.values[0]"},
                "clearable":False,
                "options":{"function":"config.adata.var.columns.values"}
            },
        },
        {
            "input":"Dropdown",
            "name":"values_to_plot",
            "description":"Values to plot.",
            "properties":{
                "value":"scores",
                "clearable":False,
                "options":['scores', 'logfoldchanges', 'pvals', 'pvals_adj']
            },
        },
        {
            "input":"Dropdown",
            "name":"plot_cluster",
            "description":"Values to plot.",
            "properties":{
                "value":{"function":"config.adata.uns['rank_genes_groups']['scores'].dtype.names[0]"},
                "clearable":False,
                "options":{"function":"config.adata.uns['rank_genes_groups']['scores'].dtype.names"}
            },
        },
    ]
}

def rank_genes_groups_f(adata, kwargs):

    sc.tl.rank_genes_groups(
                config.adata,
                groupby=kwargs["input"],
                use_raw=kwargs["use_raw"],
                layer=kwargs["layer"],
                n_genes=kwargs["n_genes"],
                method=kwargs["method"],
                corr_method=kwargs["corr_method"],
                tie_correct=kwargs["tie_correct"],
                pts=kwargs["pts"],
                )

def rank_genes_groups_plot():
    
    plot_params = get_node(config.selected)["data"]["plot"]

    data_array, labels_x, labels_y = de2array(config.adata.uns["rank_genes_groups"], plot_params["values_to_plot"], plot_params["plot_n_genes"])

    d = config.adata.var[plot_params["var_key"]]
    labels_x = d.loc[labels_x]

    if plot_params["plot_style"] == "heatmap":

        fig = plot_clustermap(data_array, labels_x, labels_y, style="heatmap")

        fig = plot_center(fig)

    elif plot_params["plot_style"] == "scattermap":

        fig = plot_clustermap(data_array, labels_x, labels_y, style="scattermap")

        fig = plot_center(fig)

    elif plot_params["plot_style"] ==  "table":

        ll = ["scores","pvals","pvals_adj","logfoldchanges"]

        df = pd.DataFrame()
        df["names"] = config.adata.var.loc[config.adata.uns["rank_genes_groups"]["names"][plot_params["plot_cluster"]],plot_params["var_key"]].values
        for i in ll:
            df[i] = config.adata.uns["rank_genes_groups"][i][plot_params["plot_cluster"]]

        fig = dag.AgGrid(
            id="table",
            rowData=df.to_dict("records"),
            columnDefs=[{"field":i} for i in df.columns.values],
            defaultColDef={"resizable": True, "sortable": True, "filter": False, "editable": False},
            columnSize="sizeToFit",
        )

    return fig

config.methods["rank_genes_groups"] = {

    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": rank_genes_groups_args,

    "function":rank_genes_groups_f,

    "plot":rank_genes_groups_plot,

}