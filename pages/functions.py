import scanpy as sc
import numpy as np
import pandas as pd

def f_qc(adata):
    if "qc_total_counts" not in adata.obs.columns.values:
        # sc.pp.calculate_qc_metrics(adata,percent_top=(3,),inplace=True)#np.round(np.int,np.linspace(1,adata.shape[1],5)))
        adata.obs["qc_total_counts"] = np.array(adata.layers["X_raw"].sum(axis=1)).reshape(-1)
        adata.obs["qc_n_genes_by_counts"] = np.array((adata.layers["X_raw"] > 0).sum(axis=1)).reshape(-1)
    if "qc" not in adata.uns.keys():
        adata.uns["qc"] = {
            "qc_total_counts":{"Minimum threshold":adata.obs["qc_total_counts"].min(), "Maximum threshold":adata.obs["qc_total_counts"].max()},
            "qc_n_genes_by_counts":{"Minimum threshold":adata.obs["qc_n_genes_by_counts"].min(), "Maximum threshold":adata.obs["qc_n_genes_by_counts"].max()},
        }

def f_qc_pattern(adata, pattern):
    if len(pattern) > 0:
        p = [i.startswith(pattern) for i in adata.var["Gene"].values]
        adata.obs["qc_"+pattern] = np.array(adata.layers["X_raw"][:,p].sum(axis=1)).reshape(-1)/adata.obs["qc_total_counts"]
        adata.uns["qc"]["qc_"+pattern] = [adata.obs["qc_"+pattern].min(), adata.obs["qc_"+pattern].max()]

def f_qc_summary(adata):
    statistics = {}
    statistics["Type"] = np.array(adata.uns["qc"].keys())
    statistics["Threshold min"] = [adata.uns["qc"][i]["Minimum threshold"] for i in adata.uns["qc"]]
    statistics["Threshold max"] = [adata.uns["qc"][i]["Maximum threshold"] for i in adata.uns["qc"]]
    statistics["Cells rejected"] = np.array([np.sum((adata.obs[i].values<adata.uns["qc"][i]["Minimum threshold"]) + (adata.obs[i].values>adata.uns["qc"][i]["Maximum threshold"]) > 0) for i in adata.uns["qc"].keys()])
    statistics["Cells rejected %"] = statistics["Cells rejected"]/adata.shape[0]

    return pd.DataFrame(statistics)