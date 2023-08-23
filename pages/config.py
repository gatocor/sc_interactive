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

graph = [
    {
        'data': {
            'id':'Raw', 
            'name':'Raw',
            'type':'Raw', 
            'method':'Raw', 
            'label':'Raw', 
            'color':'white',
            'computed':True,
            'opacity':1,
            'summary':'Raw', 
            'image':'../assets/Raw.png',
            'parameters':{'input':'Raw'}}, 
        'position':{'x':-200,'y':0},
    }
]
active_node_parameters = {}

functions = {}
functions_method = {}
functions_method_rm = {}
functions_method_rename = {}
functions_plot = {}
arguments = {}
tclick = 0
selected = 'Raw'
