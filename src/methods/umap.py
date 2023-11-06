import scanpy as sc

from general import *

umap_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"Input",
            "name":"min_dist",
            "description":"The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out. The default of in the umap-learn package is 0.1.",
            "properties" : {
                "value":.5,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"spread",
            "description":"The effective scale of embedded points. In combination with min_dist this determines how clustered/clumped the embedded points are.",
            "properties" : {
                "value":1.0,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"n_dimensions",
            "description":"The number of dimensions of the embedding.",
            "properties" : {
                "value":2,
                "clearable":False,
                "options":[2,3]
            }
        },
        {
            "input":"Input",
            "name":"maxiter",
            "description":"The number of iterations (epochs) of the optimization. Called n_epochs in the original UMAP.",
            "properties" : {
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"alpha",
            "description":"The initial learning rate for the embedding optimization.",
            "properties" : {
                "value":1.0,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"gamma",
            "description":"Weighting applied to negative samples in low dimensional embedding optimization. Values higher than one will result in greater weight being given to negative samples.",
            "properties" : {
                "value":1.0,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"negative_sample_rate",
            "description":"The number of negative edge/1-simplex samples to use per positive edge/1-simplex sample in optimizing the low dimensional embedding.",
            "properties" : {
                "value":5,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"init_pos",
            "description":"How to initialize the low dimensional embedding.",
            "properties" : {
                "value":'spectral',
                "clearable":False,
                "options":['spectral','random']
            }
        },
        {
            "input":"Input",
            "name":"random_state",
            "description":"If int, random_state is the seed used by the random number generator; If RandomState or Generator, random_state is the random number generator; If None, the random number generator is the RandomState instance used by np.random.",
            "properties" : {
                "value":0,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"a",
            "description":"More specific parameters controlling the embedding. If None these values are set automatically as determined by min_dist and spread.",
            "properties" : {
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"Input",
            "name":"b",
            "description":"More specific parameters controlling the embedding. If None these values are set automatically as determined by min_dist and spread.",
            "properties" : {
                "value":None,
                "type":"number"
            }
        },
    ],

    "postexecution" : [],

    "plot" : ARGS_COLOR

}

def umap_f(adata, kwargs):

    sc.tl.umap(adata,
            #    neighbors_key=kwargs["neighbors_key"],
                min_dist=kwargs["min_dist"],
                spread=kwargs["spread"],
                n_components=kwargs["n_dimensions"],
                maxiter=kwargs["maxiter"],
                alpha=kwargs["alpha"],
                gamma=kwargs["gamma"],
                negative_sample_rate=kwargs["negative_sample_rate"],
                init_pos=kwargs["init_pos"],
                random_state=kwargs["random_state"],
                a=kwargs["a"],
                b=kwargs["b"],
            )


def umap_plot():

    c = get_color()
    X = config.adata.obsm["X_umap"]

    n_dims = get_node(config.selected)["data"]["parameters"]["n_dimensions"]

    if n_dims == 2:

        fig = px.scatter(
                        x=X[:,0],
                        y=X[:,1],
                        color=c,
                        height=PLOTHEIGHT,
                        width=PLOTWIDTH
                )

    else:

        fig = px.scatter_3d(
                            x=X[:,0],
                            y=X[:,1],
                            z=X[:,2],
                            color=c,
                            opacity=0.8,
                            height=PLOTHEIGHT,
                            width=PLOTWIDTH
                )

    # fig.update_layout({"yaxis":{"scaleanchor":"x","scaleratio":1}})
    
    return  plot_center(dcc.Graph(figure=fig))

config.methods["umap"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":False,
    },

    "args": umap_args,

    "function": umap_f,

    "plot": umap_plot,

}