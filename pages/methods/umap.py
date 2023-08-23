import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import plotly.graph_objs as go
import dash

def umap_args(adata):

    options = []
    value = None
    if "__interactive__" in adata.uns.keys():
        options = [i for i,j in adata.uns["__interactive__"].items() if j["type"]=="Neighbors"] 
    if options != []:      
        value = options[0]

    return [
        "UMAP",
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Use the indicated neighbors.",
            "value":value,
            "clearable":False,
            "options":options
        },
        {
            "input":"Input",
            "name":"min_dist",
            "description":"The effective minimum distance between embedded points. Smaller values will result in a more clustered/clumped embedding where nearby points on the manifold are drawn closer together, while larger values will result on a more even dispersal of points. The value should be set relative to the spread value, which determines the scale at which embedded points will be spread out. The default of in the umap-learn package is 0.1.",
            "value":.5,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"spread",
            "description":"The effective scale of embedded points. In combination with min_dist this determines how clustered/clumped the embedded points are.",
            "value":1.0,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"n_components",
            "description":"The number of dimensions of the embedding.",
            "value":2,
            "clearable":False,
            "options":[2,3]
        },
        {
            "input":"Input",
            "name":"maxiter",
            "description":"The number of iterations (epochs) of the optimization. Called n_epochs in the original UMAP.",
            "value":None,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"alpha",
            "description":"The initial learning rate for the embedding optimization.",
            "value":1.0,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"gamma",
            "description":"Weighting applied to negative samples in low dimensional embedding optimization. Values higher than one will result in greater weight being given to negative samples.",
            "value":1.0,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"negative_sample_rate",
            "description":"The number of negative edge/1-simplex samples to use per positive edge/1-simplex sample in optimizing the low dimensional embedding.",
            "value":5,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"init_pos",
            "description":"How to initialize the low dimensional embedding.",
            "value":'spectral',
            "clearable":False,
            "options":['spectral','random']
        },

        {
            "input":"Input",
            "name":"random_state",
            "description":"If int, random_state is the seed used by the random number generator; If RandomState or Generator, random_state is the random number generator; If None, the random number generator is the RandomState instance used by np.random.",
            "value":0,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"a",
            "description":"More specific parameters controlling the embedding. If None these values are set automatically as determined by min_dist and spread.",
            "value":None,
            "type":"number"
        },
        {
            "input":"Input",
            "name":"b",
            "description":"More specific parameters controlling the embedding. If None these values are set automatically as determined by min_dist and spread.",
            "value":None,
            "type":"number"
        },
    ]

def f_umap(adata, name_analysis, **kwargs):

    sc.tl.umap(adata,
                neighbors_key=kwargs["input"],
                min_dist=kwargs["min_dist"],
                spread=kwargs["spread"],
                n_components=kwargs["n_components"],
                maxiter=kwargs["maxiter"],
                alpha=kwargs["alpha"],
                gamma=kwargs["gamma"],
                negative_sample_rate=kwargs["negative_sample_rate"],
                init_pos=kwargs["init_pos"],
                random_state=kwargs["random_state"],
                a=kwargs["a"],
                b=kwargs["b"],
                )

    adata.obsm["X_"+name_analysis] = adata.obsm["X_umap"]
    del adata.obsm["X_umap"]

def make_umap_plots1(adata, name_analysis):

    if "X_"+name_analysis in adata.obsm.keys() and "n_components" in adata.uns["__interactive__"][name_analysis]["params"].keys():

        if adata.uns["__interactive__"][name_analysis]["params"]["n_components"] == 2:
            x = adata.obsm["X_"+name_analysis][:,0]
            y = adata.obsm["X_"+name_analysis][:,1]

            return [
                dbc.Col(
                    dcc.Graph(
                        figure = {'data':[
                                    go.Scattergl(
                                        x=x,
                                        y=y,
                                        mode='markers',
                                    )],
                                    'layout':{
                                            'yaxis': {
                                                'scaleanchor': 'x',
                                                'scaleratio': 1
                                            },
                                            'width':900,
                                            'height':800,
                                    },
                        }
                    )
                ),
            ]
        else:
            x = adata.obsm["X_"+name_analysis][:,0]
            y = adata.obsm["X_"+name_analysis][:,1]
            z = adata.obsm["X_"+name_analysis][:,2]

            return [
                dbc.Col(
                    dcc.Graph(
                        figure = {'data':[
                                    go.Scatter3d(
                                        x=x,
                                        y=y,
                                        z=z,
                                        mode='markers',
                                        name='Min threshold',
                                    )],
                                    'layout':{
                                            'yaxis': {
                                                'scaleanchor': 'x',
                                                'scaleratio': 1
                                            },
                                            'zaxis': {
                                                'scaleanchor': 'x',
                                                'scaleratio': 1
                                            },
                                            'width':900,
                                            'height':800,
                                    },
                        }
                    )
                ),
            ]
    else:
        return []

def make_umap_plots2(adata, name_analysis):

    return []