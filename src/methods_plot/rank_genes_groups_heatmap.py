
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
import base64
from io import BytesIO

from general import *

def rank_genes_groups_heatmap_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.rank_genes_groups_heatmap(
        config.adata,
        groups=type_formater(kwargs["groups"],typing.Union[str, typing.Sequence[str]]),
        n_genes=type_formater(kwargs["n_genes"],typing.Optional[int]),
        groupby=type_formater(kwargs["groupby"],typing.Optional[str]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        var_names=type_formater(kwargs["var_names"],typing.Union[typing.Sequence[str], typing.Mapping[str, typing.Sequence[str]], str]),
        min_logfoldchange=type_formater(kwargs["min_logfoldchange"],typing.Optional[float]),
        key=type_formater(kwargs["key"],str),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Optional[bool]),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["rank_genes_groups_heatmap"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='groups', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groups'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_genes', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['n_genes'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='groupby', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groupby'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='gene_symbols', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['gene_symbols'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='var_names', 
        description="typing.Union[typing.Sequence[str], typing.Mapping[str, typing.Sequence[str]], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['var_names'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='min_logfoldchange', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['min_logfoldchange'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='key', 
        description="<class 'str'>", 
        visible=dict(function="str(None)!=config.active_plot_parameters['key'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='show', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['show'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='save', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['save'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = rank_genes_groups_heatmap_plot
)