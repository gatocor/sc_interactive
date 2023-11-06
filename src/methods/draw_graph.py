import scanpy as sc

from general import *

draw_graph_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"Dropdown",
            "name":"layout",
            "description":"'fa' (ForceAtlas2) or any valid igraph layout. Of particular interest are 'fr' (Fruchterman Reingold), 'grid_fr' (Grid Fruchterman Reingold, faster than 'fr'), 'kk' (Kamadi Kawai', slower than 'fr'), 'lgl' (Large Graph, very fast), 'drl' (Distributed Recursive Layout, pretty fast) and 'rt' (Reingold Tilford tree layout).",
            "properties" : {
                "value":'fa',
                "clearable":False,
                "options":['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa']
            }
        },
        {
            "input":"Dropdown",
            "name":"root",
            "description":"Root for tree layouts.",
            "properties" : {
                "value":None,
                "clearable":True,
                "options":{"function":"list(np.arange(config.adata.shape[0]))"}
            }
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"For layouts with random initialization like 'fr', change this to use different intial states for the optimization. If None, no seed is set.",
            "properties" : {
                "value":0,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"init_pos",
            "description":"'paga'/True, None/False, or any valid 2d-.obsm key. Use precomputed coordinates for initialization. If False/None (the default), initialize randomly.",
            "properties" : {
                "value":None,
                "clearable":True,
                "options":{"function":"['paga']+[i for i,j in config.adata.obsm.items() if j.shape[1] == 2]"}
            }
        },
    ],

    "postexecution" : [],

    "plot" : ARGS_COLOR

}

def draw_graph_f(adata, kwargs):

    sc.tl.draw_graph(adata,
        layout = kwargs["layout"],
        init_pos = kwargs["init_pos"],
        root = kwargs["root"],
        random_state = kwargs["random_state"],
    )

def draw_graph_plot():

    layout = get_node(config.selected)["data"]["parameters"]["layout"]
    c = get_color()
    X = config.adata.obsm["X_draw_graph_"+layout]

    fig = px.scatter(
                    x=X[:,0],
                    y=X[:,1],
                    color=c,
                    height=PLOTHEIGHT,
                    width=PLOTWIDTH
            )
    
    return  plot_center(dcc.Graph(figure=fig))

config.methods["draw_graph"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": draw_graph_args,

    "function": draw_graph_f,

    "plot": draw_graph_plot,

}