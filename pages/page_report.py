import dash
from dash import html
from dash import dcc
import dash_bootstrap_components as dbc
import os
import scanpy as sc
from . import config
import pandas as pd
from .functions import *

from app import app


# change to app.layout if running as single page app instead
def layout():

    return dbc.Container()

