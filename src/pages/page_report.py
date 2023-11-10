import dash
from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
import os
import scanpy as sc
import pandas as pd
import dash_editor_components
import re

from general import *

from app import app

from PIL import Image

# load_analysis("../Raw_h5ad.sc")

# change to app.layout if running as single page app instead
def layout():

    report = ""
    for i in get_nodes():
        report += i["data"]["report"]

    return dbc.Container([
        dbc.Row(
            id = "report",
            children = markdown_to_dash(report)
        ),
    ])
