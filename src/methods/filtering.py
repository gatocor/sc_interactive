import numpy as np
import pandas as pd

from general import *

from app import app

args = {
    "execution" : [
            ARGINPUT,
            {
                "input":"AgTable",
                "name":"filtering_thresholds",
                "description":"Metrics of quality control to compute",
                "properties":{
                    "header":[
                        { "headerName": ".obs", "field":"obs", "editable": True,
                        "cellEditor": "agSelectCellEditor",
                        "cellEditorParams": {"values": {"function":"filtering_names()"}},
                        },
                    ],
                    "data":{"function":"filtering_dict()"},
                },
                "addRows":{"name":""},
                "deleteRows": True,
                "recomputeAfter": ["input"] 
            },
    ],

    "postexecution" : [],

    "plot" : []
}

def filtering_names():

    return [i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in [bool]]

def filtering_dict():

    return [{"obs":i} for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in [bool] and i.endswith("--keep")]

def filtering_f(adata, kwargs):

    d = pd.DataFrame()
    keep = np.ones(adata.X.shape[0], bool)
    for i in kwargs["filtering_thresholds"]:
        name = i["obs"]
        d[name] = [np.mean(adata.obs[name].values)]
        keep = adata.obs[name].values * keep
    d["TOTAL"] = [np.mean(keep)]

    if config.selected in config.adata.uns.keys():
        config.adata.uns[config.selected]["statistics"] = d
    else:
        config.adata.uns[config.selected] = {"statistics": d}

    #remove
    config.adata = config.adata[keep,:]

    return

def filtering_plot():

    return [
        dag.AgGrid(
            id="filtering_ag_table",
            columnDefs=[
                {"headerName":i, "field":i, "valueFormatter": {"function": "`${params.value * 100}%`"}} for i in config.adata.uns[config.selected]["statistics"].columns.values
            ],
            rowData=config.adata.uns[config.selected]["statistics"].to_dict("records"),
            columnSize="sizeToFit",
        )
    ]

config.methods["filtering"] = {
    
    "properties":{
        "type":"QC",
        "make_new_h5ad":True,
    },

    "args": args,

    "function": filtering_f,

    "plot": filtering_plot,

}