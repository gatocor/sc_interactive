import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import dash
import dash_ag_grid as dag

from general import *

def args_differential_expression():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_layer = [i for i in config.adata.layers.keys()] 

    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Use the indicated neighbors.",
            "value":None,
            "clearable":False,
            "options":options
        },
        {
            "input":"BooleanSwitch",
            "name":"use_raw",
            "description":"Use raw attribute of adata if present.",
            "value":True,
        },
        {
            "input":"Dropdown",
            "name":"layer",
            "description":"Key from adata.layers whose value will be used to perform tests on.",
            "value":None,
            "clearable":True,
            "options":options_layer,
            "summary":True
        },
        {
            "input":"Input",
            "name":"n_genes",
            "description":"The number of genes that appear in the returned tables. Defaults to all genes.",
            "value":None,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"method",
            "description":"The default method is 't-test', 't-test_overestim_var' overestimates variance of each group, 'wilcoxon' uses Wilcoxon rank-sum, 'logreg' uses logistic regression. See [Ntranos18], here and here, for why this is meaningful.",
            "value":"t-test",
            "clearable":False,
            "options":["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
            "summary":True
        },
        {
            "input":"Dropdown",
            "name":"corr_method",
            "description":"p-value correction method. Used only for 't-test', 't-test_overestim_var', and 'wilcoxon'.",
            "value":"benjamini-hochberg",
            "clearable":False,
            "options":["benjamini-hochberg", "bonferroni"],
            "summary":True
        },
        {
            "input":"BooleanSwitch",
            "name":"tie_correct",
            "description":"Use tie correction for 'wilcoxon' scores. Used only for 'wilcoxon'.",
            "value":False,
        },
        {
            "input":"BooleanSwitch",
            "name":"pts",
            "description":"Compute the fraction of cells expressing the genes.",
            "value":False,
        },
    ]

def f_differential_expression(name_analysis, kwargs):

    if "log1p" not in config.adata.uns.keys(): #Solve in case there is a problem with it
        config.adata.uns["log1p"] = {"base":10}
    elif "base" not in config.adata.uns["log1p"].keys():
        config.adata.uns["log1p"] = {"base":10}

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
                key_added=name_analysis
                )
    
    pos = get_node_pos(name_analysis)
    config.graph[pos]["data"]["plotting"] = {"n_genes":2,"style":"table","cluster":config.adata.uns[config.selected]["scores"].dtype.names[0]}

def rm_differential_expression(name_analysis):

    del config.adata.uns[name_analysis]

    return

def rename_differential_expression(name_analysis, name_new_analysis):

    config.adata.uns[name_new_analysis] = config.adata.uns[name_analysis]
    del config.adata.uns[name_analysis]

    return

def plot_differential_expression(name_analysis):
    
    node = get_node(name_analysis)

    if not node["data"]["computed"]:
        return []

    l = []

    data_array, labels_x, labels_y = de2array(config.adata.uns[name_analysis], node["data"]["plotting"]["n_genes"])

    if node["data"]["plotting"]["style"] == "heatmap":
        l = [
                dcc.Dropdown(
                    id = "differential_expression_plot_style",
                    value=get_node(config.selected)["data"]["plotting"]["style"],
                    options=["heatmap","scattermap","table"],
                    clearable=False
                ),
                dcc.Dropdown(
                    id = "differential_expression_plot_n_genes",
                    value=get_node(config.selected)["data"]["plotting"]["n_genes"],
                    options=[i for i in range(1,10)],
                    clearable=False
                )
        ] + plot_clustermap(data_array, labels_x, labels_y, style="heatmap")
    elif node["data"]["plotting"]["style"] == "scattermap":
        l = [
                dcc.Dropdown(
                    id = "differential_expression_plot_style",
                    value=get_node(config.selected)["data"]["plotting"]["style"],
                    options=["heatmap","scattermap","table"],
                    clearable=False
                ),
                dcc.Dropdown(
                    id = "differential_expression_plot_n_genes",
                    value=get_node(config.selected)["data"]["plotting"]["n_genes"],
                    options=[i for i in range(1,10)],
                    clearable=False
                )
        ] + plot_clustermap(data_array, labels_x, labels_y, style="scattermap")
    elif node["data"]["plotting"]["style"] == "table":
        ll = ["names","scores","pvals","pvals_adj"]

        df = pd.DataFrame(columns = ll,
                        data={i:config.adata.uns[config.selected][i][get_node(config.selected)["data"]["plotting"]["cluster"]] for i in ll}
        )

        l = [
            dcc.Dropdown(
                id = "differential_expression_plot_style",
                value=get_node(config.selected)["data"]["plotting"]["style"],
                options=["heatmap","clustermap","table"],
                clearable=False
            ),
            dcc.Dropdown(
                id = "differential_expression_cluster",
                value=get_node(config.selected)["data"]["plotting"]["cluster"],
                options=config.adata.uns[config.selected]["scores"].dtype.names,
                clearable=False
            ),
            plot_table(df)
        ]

    return l

# @app.callback(
#     dash.Output("analysis_use_raw","on", allow_duplicate=True),
#     dash.Input("analysis_input","value"),
#     prevent_initial_call=True
# )
# def change_raw(v):

#     if config.adata.raw:
#         return True
#     else:
#         return False

# @app.callback(
#     dash.Output("analysis_plot","children", allow_duplicate=True),
#     dash.Input("differential_expression_plot_style", "value"),
#     prevent_initial_call = True
# )
# def change_style(val):

#     prevent_race("differential_expression")

#     pos = get_node_pos(config.selected)
#     config.graph[pos]["data"]["plotting"]["style"] = val

#     return plot_differential_expression(config.selected)

# @app.callback(
#     dash.Output("analysis_plot","children", allow_duplicate=True),
#     dash.Input("differential_expression_plot_n_genes", "value"),
#     prevent_initial_call = True
# )
# def change_style(val):

#     prevent_race("differential_expression")

#     pos = get_node_pos(config.selected)
#     config.graph[pos]["data"]["plotting"]["n_genes"] = val

#     return plot_differential_expression(config.selected)

# @app.callback(
#     dash.Output("analysis_plot","children", allow_duplicate=True),
#     dash.Input("differential_expression_cluster", "value"),
#     prevent_initial_call = True
# )
# def change_style(val):

#     prevent_race("differential_expression")

#     pos = get_node_pos(config.selected)
#     config.graph[pos]["data"]["plotting"]["cluster"] = val

#     return plot_differential_expression(config.selected)