import scanpy as sc
import pandas as pd
import os
from .functions import *

CACHEFOLDER="__h5adcache__/"

def make_file_path(name):
    return "./__graph__/"+name.split("/")[-1].split(".")[0]+"_qc"+".h5ad"

# Initialize the global variable
folder_path = './'  # Replace with the actual path to the folder containing .h5ad files
h5ad_files = []
for folder, dirs, files in os.walk("."):
    for file in files:
        if ".h5ad" in file:
            h5ad_files.append(folder+"/"+file)
selected_file = None
file_path = h5ad_files[0] #None
try: 
    f = make_file_path(file_path)
    adata = sc.read(f) #None
except:
    adata = sc.read(file_path) #None

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
            'parameters':{'input':None}}, 
        'position':{'x':-200,'y':0},
    }
]
active_node_parameters = {}
max_x = 0

functions = {}
functions_method = {}
functions_method_rm = {}
functions_method_rename = {}
functions_plot = {}
arguments = {}
tclick = 0
selected = 'Raw'
