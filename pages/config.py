import scanpy as sc
import pandas as pd
import os
from .functions import *

# Initialize the global variable
folder_path = './'  # Replace with the actual path to the folder containing .h5ad files
h5ad_files = np.sort([f for f in os.listdir(folder_path) if f.endswith('.h5ad')])
old_selected_file = h5ad_files[0]
file_path = os.path.join(folder_path, h5ad_files[0])
adata = sc.read(file_path)

#QC
# Perform necessary computations or filtering on the adata object
f_qc(adata)

#Options
qc_control_options = [i for i in adata.obs.columns.values if i.startswith("qc_")]

#Statistics
qc_df_statistics = pd.DataFrame(adata.uns["qc"])

#Patterns
qc_dict_patterns = {
    'Concept': [],
    'Pattern': [],
    'Genes': []
}
qc_df_patterns = pd.DataFrame(qc_dict_patterns)

qc_dict_threshold = {
    'Control measure': [], 
    'Minimum Threshold': [], 
    'Maximum Threshold': []
}
qc_df_threshold = pd.DataFrame(qc_dict_threshold)

qc_n_clicks_old = 0

qc_summary = f_qc_summary(adata)
