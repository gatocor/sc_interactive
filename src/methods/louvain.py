import scanpy as sc

from general import *

louvain_args = {

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
            "input":"Dropdown",
            "name":"flavor",
            "description":"Choose between to packages for computing the clustering.",
            "properties" : {
                "value":"vtraag",
                "clearable":False,
                "options":["vtraag","igraph","rapids"],
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
    ],

    "postexecution" : [],

    "plot" : ARGS_REPRESENTATION
}

def louvain_f(adata, kwargs):

    sc.tl.louvain(
        adata,
        resolution=kwargs["resolution"],
        flavor=kwargs["flavor"],
        random_state=kwargs["random_state"], 
        directed=kwargs["directed"], 
        use_weights=kwargs["use_weights"], 
    )

def louvain_plot():

    c = config.adata.obs["louvain"]

    fig = get_representation(color=c)
    fig.update_layout(height=1200, width=1200, autosize=True, showlegend=False)

    return plot_center(dcc.Graph(figure=fig))

config.methods["louvain"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": louvain_args,

    "function": louvain_f,

    "plot": louvain_plot,

}
