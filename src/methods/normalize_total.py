import scanpy as sc

from general import *

normalize_total_args = {

    "execution" : [
        ARGINPUT,
        {
            "input":"Input",
            "name":"target_sum",
            "description":"If None, after normalization, each observation (cell) has a total count equal to the median of total counts for observations (cells) before normalization.",
            "properties" : {
                "value":None,
                "type":"number"
            }
        },
        {
            "input":"BooleanSwitch",
            "name":"exclude_highly_expressed",
            "description":"Exclude (very) highly expressed genes for the computation of the normalization factor (size factor) for each cell. A gene is considered highly expressed, if it has more than max_fraction of the total counts in at least one cell. The not-excluded genes will sum up to target_sum.",
            "properties" : {
                "on":False,
            }
        },
        {
            "input":"Input",
            "name":"max_fraction",
            "description":"If exclude_highly_expressed=True, consider cells as highly expressed that have more counts than max_fraction of the original total counts in at least one cell.",
            "properties":{
                "value":0.05,
                "type":"number"
            },
            "visible" : {"function":"config.active_node_parameters['exclude_highly_expressed']"}
        },
    ],

    "postexecution" : [],

    "plot" : [],

}

def normalize_total_f(adata, kwargs):
        
    sc.pp.normalize_total(config.adata,
        target_sum=kwargs["target_sum"],
        exclude_highly_expressed=kwargs["exclude_highly_expressed"],
        max_fraction=kwargs["max_fraction"]
    )

def normalize_total_plot():

    return []

config.methods["normalize_total"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":True,
    },

    "args": normalize_total_args,

    "function": normalize_total_f,

    "plot": normalize_total_plot,

}