
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

def rank_genes_groups_violin_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.rank_genes_groups_violin(
        config.adata,
        groups=type_formater(kwargs["groups"],typing.Optional[typing.Sequence[str]]),
        n_genes=type_formater(kwargs["n_genes"],int),
        gene_names=type_formater(kwargs["gene_names"],typing.Optional[typing.Iterable[str]]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        key=type_formater(kwargs["key"],typing.Optional[str]),
        split=type_formater(kwargs["split"],bool),
        scale=type_formater(kwargs["scale"],str),
        strip=type_formater(kwargs["strip"],bool),
        jitter=type_formater(kwargs["jitter"],typing.Union[int, float, bool]),
        size=type_formater(kwargs["size"],int),
        ax=ax,
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["rank_genes_groups_violin"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='groups', 
        description="typing.Optional[typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groups'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_genes', 
        description="<class 'int'>", 
        visible=dict(function="str(20)!=config.active_plot_parameters['n_genes'] or config.show_plot"),
        properties=dict(value="20",type="text")
    ),
    dict(
        input='Input', 
        name='gene_names', 
        description="typing.Optional[typing.Iterable[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['gene_names'] or config.show_plot"),
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
        name='use_raw', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['use_raw'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['key'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='split', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['split'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='scale', 
        description="<class 'str'>", 
        visible=dict(function="'width'!=eval(config.active_plot_parameters['scale']) or config.show_plot"),
        properties=dict(value="'width'",type="text")
    ),
    dict(
        input='Input', 
        name='strip', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['strip'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='jitter', 
        description="typing.Union[int, float, bool]", 
        visible=dict(function="str(True)!=config.active_plot_parameters['jitter'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='size', 
        description="<class 'int'>", 
        visible=dict(function="str(1)!=config.active_plot_parameters['size'] or config.show_plot"),
        properties=dict(value="1",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = rank_genes_groups_violin_plot
)
