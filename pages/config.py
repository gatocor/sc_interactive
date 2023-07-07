import scanpy as sc
import pandas as pd
import os

# Initialize the global variable
folder_path = './'  # Replace with the actual path to the folder containing .h5ad files
h5ad_files = [f for f in os.listdir(folder_path) if f.endswith('.h5ad')]
old_selected_file = h5ad_files[0]
file_path = os.path.join(folder_path, h5ad_files[0])
adata = sc.read(file_path)
# Perform necessary computations or filtering on the adata object
if "total_counts" not in adata.obs.columns.values:
    sc.pp.calculate_qc_metrics(adata,percent_top=(3,),inplace=True)#np.round(np.int,np.linspace(1,adata.shape[1],5)))
if "qc" not in adata.uns.keys():
    adata.uns["qc"] = {
        "total_counts":[adata.obs["total_counts"].min(), adata.obs["total_counts"].max()],
        "n_genes_by_counts":[adata.obs["n_genes_by_counts"].min(), adata.obs["n_genes_by_counts"].max()],
        # "mt_fraction":[adata.obs["mt_fraction"].min(), adata.obs["mt_fraction"].max()]
    }

df_qc = pd.DataFrame(adata.uns["qc"])