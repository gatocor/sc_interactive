import scanpy as sc

from general import *

leiden_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"Input",
            "name":"resolution",
            "description":"A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters. Set to None if overriding partition_type to one that doesn’t accept a resolution_parameter.",
            "properties" : {
                "value":1,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"Change the initialization of the optimization.",
            "properties" : {
                "value":0,
                "type":"number"
            }
        },
        {
            "input":"BooleanSwitch",
            "name":"directed",
            "description":"Whether to treat the graph as directed or undirected.",
            "properties" : {
                "on":True,
            }
        },
        {
            "input":"BooleanSwitch",
            "name":"use_weights",
            "description":"If True, edge weights from the graph are used in the computation (placing more emphasis on stronger edges).",
            "properties" : {
                "on":True,
            }
        },
        {
            "input":"Input",
            "name":"n_iterations",
            "description":"How many iterations of the Leiden clustering algorithm to perform. Positive values above 2 define the total number of iterations to perform, -1 has the algorithm run until it reaches its optimal clustering.",
            "properties" : {
                "value":-1,
                "type":"number"
            }
        },
    ],

    "postexecution" : [],

    "plot" : ARGS_REPRESENTATION
}

def leiden_f(adata, kwargs):

    sc.tl.leiden(
        adata,
        resolution=kwargs["resolution"],
        random_state=kwargs["random_state"], 
        directed=kwargs["directed"], 
        use_weights=kwargs["use_weights"], 
        n_iterations=kwargs["n_iterations"], 
    )


def leiden_plot():

    c = config.adata.obs["leiden"]

    fig = get_representation(color=c)
    fig.update_layout(height=1200, width=1200, autosize=True, showlegend=False)

    return plot_center(dcc.Graph(figure=fig))

config.methods["leiden"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": leiden_args,

    "function": leiden_f,

    "plot": leiden_plot,

}

