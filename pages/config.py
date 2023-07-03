import scanpy as sc
import os

# Initialize the global variable
folder_path = './'  # Replace with the actual path to the folder containing .h5ad files
h5ad_files = [f for f in os.listdir(folder_path) if f.endswith('.h5ad')]
old_selected_file = h5ad_files[0]
file_path = os.path.join(folder_path, h5ad_files[0])
adata = sc.read(file_path)