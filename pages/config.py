import scanpy as sc
import pandas as pd
import os
from .functions import *

CACHEFOLDER="__h5adcache__/"

# Initialize the global variable
folder_path = './'  # Replace with the actual path to the folder containing .h5ad files
h5ad_files = np.sort([f for f in os.listdir(folder_path) if f.endswith('.h5ad')])
selected_file = None
file_path = os.path.join(folder_path, h5ad_files[0]) #None
adata = sc.read(os.path.join(folder_path, h5ad_files[0])) #None
cache = False

max_scatter_points = 10000

#Options
qc_summary_x = "total_counts"
qc_summary_y = "n_genes_by_counts"
qc_summary_color = "total_counts"

#Statistics
#qc_df_statistics = pd.DataFrame(adata.uns["qc"])

#Patterns
qc_dict_patterns = {
    'Concept': [],
    'Pattern': [],
    'Genes': []
}
qc_df_patterns = pd.DataFrame(qc_dict_patterns)

qc_dict_threshold = {
    'Control measure': [], 
    'Minimum threshold': [], 
    'Maximum threshold': []
}
qc_df_threshold = pd.DataFrame(qc_dict_threshold)

qc_summary = f_qc_summary(adata)
