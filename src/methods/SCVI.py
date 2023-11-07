import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scvi
scvi.settings.seed = 0

from general import *

#Only show color when plot type is components
ARGS_COLOR_scvi = deepcopy(ARGS_COLOR)
ARGS_COLOR_scvi[0]["visible"] = {"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['components']"}

scvi_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"Dropdown",
            "name":"layer",
            "description":"(Optional[str] (default: None)) - if not None, uses this as the key in adata.layers for raw count data.",
            "properties":{
                "value":None,
                "clearable":True,
                "options":{"function":"['X']+[i for i in config.adata.layers.keys()]"}
            }
        },
        {
            "input":"Dropdown",
            "name":"batch_key",
            "description":"(Optional[str] (default: None)) - key in adata.obs for batch information. Categories will automatically be converted into integer categories and saved to adata.obs['_scvi_batch']. If None, assigns the same batch to all the data.",
            "properties":{
                "value":None,
                "clearable":True,
                "options":{"function":"[i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in ['str','int','categorical']]"}
            }
        },
        {
            "input":"Dropdown",
            "name":"labels_key",
            "description":"(Optional[str] (default: None)) - key in adata.obs for label information. Categories will automatically be converted into integer categories and saved to adata.obs['_scvi_labels']. If None, assigns the same label to all the data.",
            "properties":{
                "value":None,
                "clearable":True,
                "options":{"function":"[i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in ['str','int','categorical']]"}
            }
        },
        {
            "input":"Dropdown",
            "name":"size_factor_key",
            "description":"(Optional[str] (default: None)) - key in adata.obs for size factor information. Instead of using library size as a size factor, the provided size factor column will be used as offset in the mean of the likelihood. Assumed to be on linear scale.",
            "properties":{
                "value":None,
                "clearable":True,
                "options":{"function":"[i for i in config.adata.obs.columns.values]"}
            }
        },
        {
            "input":"AgTable",
            "name":"categorical_covariate_keys",
            "description":"(Optional[List[str]] (default: None)) - keys in adata.obs that correspond to categorical data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.",
            "properties":{
                "header":[
                    { "headerName": ".obs", "field":"obs", "editable": True,
                    "cellEditor": "agSelectCellEditor",
                    "cellEditorParams": {"values": {"function":"list(config.adata.obs.columns.values)"}},
                    },
                ],
                "data":[
                ],
            },
            "addRows":{"obs":""},
            "deleteRows": True                
        },
        {
            "input":"AgTable",
            "name":"continuous_covariate_keys",
            "description":"(Optional[List[str]] (default: None)) - keys in adata.obs that correspond to continuous data. These covariates can be added in addition to the batch covariate and are also treated as nuisance factors (i.e., the model tries to minimize their effects on the latent space). Thus, these should not be used for biologically-relevant factors that you do _not_ want to correct for.",
            "properties":{
                "header":[
                    { "headerName": ".obs", "field":"obs", "editable": True,
                    "cellEditor": "agSelectCellEditor",
                    "cellEditorParams": {"values": {"function":"list(config.adata.obs.columns.values)"}},
                    },
                ],
                "data":[
                ],
            },
            "addRows":{"obs":""},
            "deleteRows": True                
        },
        {
            "input":"Dropdown",
            "name":"n_hidden",
            "description":"(int (default: 128)) - Number of nodes per hidden layer.",
            "properties":{
                "value":128,
                "clearable":False,
                "options":[2**i for i in range(10)]
            }
        },
        {
            "input":"Dropdown",
            "name":"n_latent",
            "description":"(int (default: 10)) - Dimensionality of the latent space.",
            "properties":{
                "value":10,
                "clearable":False,
                "options":[i for i in range(100)]
            }
        },
        {
            "input":"Dropdown",
            "name":"n_layers",
            "description":"(int (default: 1)) - Number of hidden layers used for encoder and decoder NNs.",
            "properties":{
                "value":1,
                "clearable":False,
                "options":[i for i in range(100)]
            }
        },
        {
            "input":"Input",
            "name":"dropout_rate",
            "description":"(float (default: 0.1)) - Dropout rate for neural networks.",
            "properties":{
                "value":.1,
                "type":"number"
            }
        },
        {
            "input":"Dropdown",
            "name":"dispersion",
            "description":"'gene' - dispersion parameter of NB is constant per gene across cells. 'gene-batch' - dispersion can differ between different batches. 'gene-label' - dispersion can differ between different labels. 'gene-cell' - dispersion can differ for every gene in every cell",
            "properties":{
                "value":"gene",
                "clearable":False,
                "options":['gene', 'gene-batch', 'gene-label', 'gene-cell']
            }
        },
        {
            "input":"Dropdown",
            "name":"gene_likelihood",
            "description":"'nb' - Negative binomial distribution. 'zinb' - Zero-inflated negative binomial distribution. 'poisson' - Poisson distribution",
            "properties":{
                "value":"zinb",
                "clearable":False,
                "options":['zinb', 'nb', 'poisson']
            }
        },
        {
            "input":"Dropdown",
            "name":"latent_distribution",
            "description":"'normal' - Normal distribution. 'ln' - Logistic normal distribution (Normal(0, I) transformed by softmax)",
            "properties":{
                "value":"normal",
                "clearable":False,
                "options":['normal', 'ln']
            }
        },
        {
            "input":"Input",
            "name":"max_epochs",
            "description":"(Optional[int] (default: None)) - Number of passes through the dataset. If None, defaults to np.min([round((20000 / n_cells) * 400), 400])",
            "properties":{
                "value":None,
                "type":"number"
            },
        },
        {
            "input":"Dropdown",
            "name":"accelerator",
            "description":"(str (default: 'auto')) - Supports passing different accelerator types ('cpu', 'gpu', 'tpu', 'ipu', 'hpu', 'mp', 'auto') as well as custom accelerator instances.",
            "properties":{
                "value":"auto",
                "options":['cpu', 'gpu', 'tpu', 'ipu', 'hpu', 'mp', 'auto']
            },
        },
        {
            "input":"Input",
            "name":"devices",
            "description":"(Union[int, List[int], str] (default: 'auto')) - The devices to use. Can be set to a non-negative index (int or str), a sequence of device indices (list or comma-separated str), the value -1 to indicate all available devices, or 'auto' for automatic selection based on the chosen accelerator. If set to 'auto' and accelerator is not determined to be 'cpu', then devices will be set to the first available device.",
            "properties":{
                "value":"auto"
            },
        },
        {
            "input":"Input",
            "name":"train_size",
            "description":"(float (default: 0.9)) - Size of training set in the range [0.0, 1.0].",
            "properties":{
                "value":0.9,
                "type":"number"
            },
        },
        {
            "input":"Input",
            "name":"validation_size",
            "description":"(Optional[float] (default: None)) - Size of the test set. If None, defaults to 1 - train_size. If train_size + validation_size < 1, the remaining cells belong to a test set.",
            "properties":{
                "value":.1,
                "type":"number"
            },
        },
        {
            "input":"BooleanSwitch",
            "name":"shuffle_set_split",
            "description":"(bool (default: True)) - Whether to shuffle indices before splitting. If False, the val, train, and test set are split in the sequential order of the data according to validation_size and train_size percentages.",
            "properties":{
                "on":True
            },
        },
        {
            "input":"Dropdown",
            "name":"batch_size",
            "description":"(int (default: 128)) - Minibatch size to use during training.",
            "properties":{
                "value":128,
                "options":[2**i for i in range(10)]
            },
        },
        {
            "input":"BooleanSwitch",
            "name":"early_stopping",
            "description":"(bool (default: False)) - Perform early stopping. Additional arguments can be passed in **kwargs. See Trainer for further options.",
            "properties":{
                "on":False
            },
        }
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
                "options":["components"] 
            }
        },
        {
            "input":"Dropdown",
            "name":"n_plot_components",
            "description":"Number of components displayed.",
            "properties":{
                "value":{"function":"min(3,config.adata.obsm['X_scVI'].shape[1])"},
                "clearable":False,
                "options":{"function":"[i for i in range(2,config.adata.obsm['X_scVI'].shape[1])]"} 
            },
            "visible":{"function":"get_node(config.selected)['data']['plot']['plot_type'] in ['components','correlation_matrix']"}
        },
    ]+ARGS_COLOR_scvi
}

def scvi_f(adata, kwargs):
        
    scvi.model.SCVI.setup_anndata(config.adata, 
        layer = kwargs["layer"],
        batch_key = kwargs["batch_key"],
        labels_key = kwargs["labels_key"],
        size_factor_key = kwargs["size_factor_key"],
        categorical_covariate_keys = [i["obs"] for i in kwargs["categorical_covariate_keys"]],
        continuous_covariate_keys = [i["obs"] for i in kwargs["continuous_covariate_keys"]],
    )
    vae = scvi.model.SCVI(adata,
        n_hidden = kwargs["n_hidden"],
        n_latent = kwargs["n_latent"],
        n_layers = kwargs["n_layers"],
        dropout_rate = kwargs["dropout_rate"],
        dispersion = kwargs["dispersion"],
        gene_likelihood = kwargs["gene_likelihood"],
        latent_distribution = kwargs["latent_distribution"],
    )
    vae.train(
        max_epochs = kwargs["max_epochs"],
        accelerator = kwargs["accelerator"],
        devices = kwargs["devices"],
        train_size = kwargs["train_size"],
        validation_size = kwargs["validation_size"],
        shuffle_set_split = kwargs["shuffle_set_split"],
        batch_size = kwargs["batch_size"],
        early_stopping = kwargs["early_stopping"],
    )
    config.adata.obsm["X_scVI"] = vae.get_latent_representation()

def scvi_plot():

    plot_params = get_node(config.selected)['data']['plot']
    plot_type = plot_params['plot_type']

    if plot_type == "components" :

        X = config.adata.obsm["X_scVI"]

        c = get_color()

        fig = make_subplots(rows=plot_params['n_plot_components']-1, cols=plot_params['n_plot_components']-1, 
                    shared_yaxes=True, 
                    shared_xaxes=True,                     
        )
        for i in range(plot_params['n_plot_components']-1):

            for j in range(i+1,plot_params['n_plot_components']):
                x_scvi = X[:,i]
                y_scvi = X[:,j]

                fig.add_traces(
                        list(px.scatter(
                                    x=x_scvi,
                                    y=y_scvi,
                                    color=c,
                        ).select_traces()),
                        rows=j, cols=i+1
                )

        fig.update_layout(width=PLOTWIDTH, height=PLOTHEIGHT, autosize=True, showlegend=True)

        return plot_center(dcc.Graph(figure=fig))

# config.methods["scvi"] = {
    
#     "properties":{
#         "type":"QC",
#         "make_new_h5ad":False,
#     },

#     "args": scvi_args,

#     "function": scvi_f,

#     "plot": scvi_plot,

# }