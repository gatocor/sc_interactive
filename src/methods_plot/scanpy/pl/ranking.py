
import numpy
from numpy import inf
import scanpy
import pandas
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scipy
import plotly.tools as tls
import cycler
import matplotlib      # pip install matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

from general import *

def ranking_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.ranking(
        config.adata,
        attr=type_formater(kwargs["attr"],typing.Literal['var', 'obs', 'uns', 'varm', 'obsm']),
        keys=type_formater(kwargs["keys"],typing.Union[str, typing.Sequence[str]]),
    )


    return fig

config.methods_plot["ranking"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='attr', 
        description="typing.Literal['var', 'obs', 'uns', 'varm', 'obsm']", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
    ),
    dict(
        input='Input', 
        name='keys', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
    ),],

    function = ranking_plot,

    docs = sc.pl.umap.__doc__
)
