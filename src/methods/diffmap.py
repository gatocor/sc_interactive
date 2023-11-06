import scanpy as sc

from general import *

ARGS_COLOR_DIFFMAP = deepcopy(ARGS_COLOR)
ARGS_COLOR_DIFFMAP[0]["visible"] = {"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['components']"}

diffmap_args = {

    "execution" : [
        ARGINPUT,

        {
            "input":"Input",
            "name":"n_comps",
            "description":"The number of dimensions of the representation.",
            "properties" : {
                "value":15,
                "type":"number"
            },
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"A numpy random seed.",
            "properties" : {
                "value":0,
            },
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
                "options":["eigenvalues", "components"] 
            }
        },
        {
            "input":"Dropdown",
            "name":"n_plot_components",
            "description":"Number of components displayed.",
            "properties":{
                "value":{"function":"min(3,len(config.adata.uns['diffmap_evals']))"},
                "clearable":False,
                "options":{"function":"list(range(1,len(config.adata.uns['diffmap_evals'])-1))"} 
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['components']"}
        },
    ]+ARGS_COLOR_DIFFMAP,

}

def diffmap_f(adata, kwargs):

    sc.tl.diffmap(adata,
        n_comps = kwargs["n_comps"],
        random_state = kwargs["random_state"],
    )

def diffmap_plot():

    plot_params = get_node(config.selected)['data']['plot']
    plot_type = plot_params['plot_type']

    if plot_type == "eigenvalues":

        y = config.adata.uns['diffmap_evals']

        fig = px.line(x=np.arange(0,len(y)),y=y)
        fig.layout["xaxis"]["title"] = "diffmap components"
        fig.layout["yaxis"]["title"] = "diffmap egeinvalues"

        return plot_center(dcc.Graph(figure=fig))

    elif plot_type == "components" :

        X = config.adata.obsm["X_diffmap"]

        c = get_color()

        fig = make_subplots(rows=plot_params['n_plot_components']-1, cols=plot_params['n_plot_components']-1, 
                            shared_yaxes=True, 
                            shared_xaxes=True,
                )
        for i in range(1,plot_params['n_plot_components']):

            for j in range(i+1,plot_params['n_plot_components']+1):
                x_pca = X[:,i]
                y_pca = X[:,j]

                fig.add_traces(
                        list(px.scatter(
                                    x=x_pca,
                                    y=y_pca,
                                    color=c,
                        ).select_traces()),
                        rows=j-1, cols=i
                )

        fig.update_layout(height=PLOTHEIGHT, width=PLOTWIDTH, autosize=True, showlegend=False)

        return plot_center(dcc.Graph(figure=fig))

config.methods["diffmap"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": diffmap_args,

    "function": diffmap_f,

    "plot": diffmap_plot,

}