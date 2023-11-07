import numpy as np
import scanpy as sc
import plotly.graph_objs as go
from plotly.subplots import make_subplots

from general import *

#Only show color when plot type is components
ARGS_COLOR_PCA = deepcopy(ARGS_COLOR)
ARGS_COLOR_PCA[0]["visible"] = {"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['components']"}

pca_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"BooleanSwitch",
            "name":"use_highly_variable",
            "description":"Whether to use highly variable genes only.",
            "properties":{
                "on":True,
            }
        },
        {
            "input":"Input",
            "name":"n_comps",
            "description":"Number of principal components to compute. Defaults to 50, or 1 - minimum dimension size of selected representation.",
            "properties":{
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"BooleanSwitch",
            "name":"zero_center",
            "description":"If True, compute standard PCA from covariance matrix. If False, omit zero-centering variables (uses TruncatedSVD), which allows to handle sparse input efficiently. Passing None decides automatically based on sparseness of the data.",
            "properties":{
                "on":True,
            }
        },
        {
            "input":"Dropdown",
            "name":"svd_solver",
            "description":"Efficient computation of the principal components of a sparse matrix currently only works with the 'arpack' or 'lobpcg' solvers.",
            "properties":{
                "value":"arpack",
                "clearable":False,
                "options":["arpack", "randomized", "auto","lobpcg"] 
            }
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change to use different initial states for the optimization.",
            "properties":{
                "value":0,
                "type":"number"
            }
        },
    ],

    "postexecution" : [],

    "plot" : [
        {
            "input":"Dropdown",
            "name":"plot_type",
            "description":"Choose the visualization aimed.",
            "properties":{
                "value":"components",
                "clearable":False,
                "options":["variance_ratio", "components", "correlation_matrix"] 
            }
        },
        {
            "input":"Dropdown",
            "name":"n_plot_components",
            "description":"Number of components displayed.",
            "properties":{
                "value":{"function":"min(3,len(config.adata.uns['pca']['variance']))"},
                "clearable":False,
                "options":{"function":"n_pcas()"} 
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['components','correlation_matrix']"}
        },
        {
            "input":"BooleanSwitch",
            "name":"show_data",
            "description":"Show scatter data.",
            "properties":{
                "on":True,
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['plot_type'] == 'components'"}
        },
        {
            "input":"BooleanSwitch",
            "name":"show_loadings",
            "description":"Show loadings.",
            "properties":{
                "on":False,
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['plot_type'] == 'components'"}
        },
        {
            "input":"Dropdown",
            "name":"n_plot_loadings",
            "description":"How many of the most important loadings to display.",
            "properties":{
                "value":{"function":"min(10,config.adata.varm['PCs'].shape[0])"},
                "clearable":False,
                "options":{"function":"n_loadings()"} 
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['components']"}
        },
        {
            "input":"Dropdown",
            "name":"n_plot_loadings_per_pc",
            "description":"How many of the most important loadings to display per component displayed.",
            "properties":{
                "value":{"function":"min(10,config.adata.varm['PCs'].shape[0])"},
                "clearable":False,
                "options":{"function":"n_loadings()"} 
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['correlation_matrix']"}
        },
        {
            "input":"Input",
            "name":"loadings_scale",
            "description":"Scale of loading arrows.",
            "properties":{
                "value":100,
                "type":"number"
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['show_loadings']"}
        },
        {
            "input":"Dropdown",
            "name":"var_key",
            "description":".var key to use as label for the loadings.",
            "properties":{
                "value":{"function":"config.adata.var.columns.values[0]"},
                "clearable":False,
                "options":{"function":"config.adata.var.columns.values"} 
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['show_loadings']"}
        },
    ]+ARGS_COLOR_PCA
}

def n_pcas():

    return list(range(1,len(config.adata.uns['pca']['variance'])))

def n_loadings():

    return list(range(1,100))

def pca_f(adata, kwargs):
        
    sc.pp.pca(config.adata,
              use_highly_variable=kwargs["use_highly_variable"],
              n_comps=kwargs["n_comps"],
              zero_center=kwargs["zero_center"],
              svd_solver=kwargs["svd_solver"],
              random_state=kwargs["random_state"]
              )

def pca_plot():

    plot_params = get_node(config.selected)['data']['plot']
    plot_type = plot_params['plot_type']

    if plot_type == "variance_ratio":

        y = config.adata.uns['pca']['variance_ratio']

        fig = px.line(x=np.arange(1,len(y)+1),y=y)
        fig.layout["xaxis"]["title"] = "pcs"
        fig.layout["yaxis"]["title"] = "variance_ratio"

        return plot_center(dcc.Graph(figure=fig))

    elif plot_type == "correlation_matrix":

        if config.adata.uns['pca']['params']['use_highly_variable']:
            X = config.adata.varm["PCs"][config.adata.var["highly_variable"].values,:plot_params['n_plot_components']]
            genes = config.adata.var[plot_params["var_key"]].values[config.adata.var["highly_variable"].values]
        else:
            X = config.adata.varm["PCs"][:,:plot_params['n_plot_components']]
            genes = config.adata.var[plot_params["var_key"]].values

        gene_names = {}
        for i in range(plot_params['n_plot_components']):
            j = np.argsort(-np.abs(X[:,i]))[:plot_params["n_plot_loadings_per_pc"]]
            for k in j:
                gene_names[k] = str(genes[k])

        m = [i for i in gene_names.keys()]
        y = [j for i,j in gene_names.items()]
        x = [str(i) for i in np.arange(1,plot_params['n_plot_components']+1)]

        fig = px.imshow(np.transpose(X[m,:]),x=y,y=x)
        fig.layout["xaxis"]["title"] = "original variables"
        fig.layout["yaxis"]["title"] = "pcs"

        return plot_center(dcc.Graph(figure=fig))

    elif plot_type == "components" :

        X = config.adata.obsm["X_pca"]

        c = get_color()

        # Loadings
        if config.adata.uns['pca']['params']['use_highly_variable']:
            X_loadings = config.adata.varm["PCs"][config.adata.var["highly_variable"].values,:plot_params['n_plot_components']]
            genes = config.adata.var[plot_params["var_key"]].values[config.adata.var["highly_variable"].values]
        else:
            X_loadings = config.adata.varm["PCs"][:,:plot_params['n_plot_components']]
            genes = config.adata.var[plot_params["var_key"]].values

        fig = make_subplots(rows=plot_params['n_plot_components']-1, cols=plot_params['n_plot_components']-1, 
                    shared_yaxes=True, 
                    shared_xaxes=True,                     
        )
        for i in range(plot_params['n_plot_components']-1):

            for j in range(i+1,plot_params['n_plot_components']):
                x_pca = X[:,i]
                y_pca = X[:,j]

                if plot_params["show_data"]:
                    if c:
                        fig.add_traces(
                                list(px.scatter(
                                            x=x_pca,
                                            y=y_pca,
                                            color=c,
                                ).select_traces()),
                                rows=j, cols=i+1
                        )
                    else:
                        fig.add_traces(
                                list(px.scatter(
                                            x=x_pca,
                                            y=y_pca,
                                ).select_traces()),
                                rows=j, cols=i+1
                        )

                if plot_params["show_loadings"]:
                    feature = {}
                    order = np.argsort(-np.sum(np.power(X_loadings[:,[i,j]],2),axis=1))[:plot_params["n_plot_loadings"]]
                    for k in order:
                        feature[k] = str(genes[k])

                    scale = plot_params["loadings_scale"]
                    for k, feature in feature.items():
                        fig.add_trace(
                            go.Scatter(
                                    x=[0,X_loadings[k, i]*scale],
                                    y=[0,X_loadings[k, j]*scale],
                                    marker={
                                        "color":"black",
                                        "symbol": "arrow-bar-up", 
                                        "angleref":"previous"
                                    }
                            ),
                            row=j, col=i+1
                        )
                        fig.add_annotation(
                            x=X_loadings[k, i]*scale,
                            y=X_loadings[k, j]*scale,
                            ax=0, ay=0,
                            xanchor="center",
                            yanchor="bottom",
                            text=feature,
                            xshift=X_loadings[k, i]*scale,
                            yshift=X_loadings[k, j]*scale,
                            row=j, col=i+1
                        )

        fig.update_layout(width=PLOTWIDTH, height=PLOTHEIGHT, autosize=True, showlegend=False)

        return plot_center(dcc.Graph(figure=fig))

config.methods["pca"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": pca_args,

    "function": pca_f,

    "plot": pca_plot,

}