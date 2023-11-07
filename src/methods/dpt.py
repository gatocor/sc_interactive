
import numpy
from numpy import inf
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scanpy as sc
import louvain
import scipy
import leidenalg
import plotly.tools as tls
import cycler
import matplotlib      # pip install matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO
from general import *

dpt_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='n_dcs', 
        description="<class 'int'>", 
        visible=dict(function="str(10)!=get_node(config.selected)['data']['parameters']['n_dcs']"),
        properties=dict(value="10",type="text")
    ),
    dict(
        input='Input', 
        name='n_branchings', 
        description="<class 'int'>", 
        visible=dict(function="str(0)!=get_node(config.selected)['data']['parameters']['n_branchings']"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='min_group_size', 
        description="<class 'float'>", 
        visible=dict(function="str(0.01)!=get_node(config.selected)['data']['parameters']['min_group_size']"),
        properties=dict(value="0.01",type="text")
    ),
    dict(
        input='Input', 
        name='allow_kendall_tau_shift', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['allow_kendall_tau_shift']"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['neighbors_key']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['copy']"),
        properties=dict(value="False",type="text")
    ),],
    postexecution = [],
    plot = [
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[str, matplotlib.colors.Colormap, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['color_map']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='palette', 
        description="typing.Union[typing.Sequence[str], cycler.Cycler, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['palette']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='show', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['show']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='save', 
        description="typing.Union[bool, str, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['save']"),
        properties=dict(value="None",type="text")
    ),]
)

def dpt_f(adata,kwargs):

    sc.tl.dpt(
        adata,
        n_dcs=type_formater(kwargs["n_dcs"],int),
        n_branchings=type_formater(kwargs["n_branchings"],int),
        min_group_size=type_formater(kwargs["min_group_size"],float),
        allow_kendall_tau_shift=type_formater(kwargs["allow_kendall_tau_shift"],bool),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def dpt_plot():

    kwargs = get_node(config.selected)['data']['plot']
    
    fig = sc.pl.dpt_groups_pseudotime(
        config.adata,
        color_map=type_formater(kwargs["color_map"],typing.Union[str, matplotlib.colors.Colormap, str]),
        palette=type_formater(kwargs["palette"],typing.Union[typing.Sequence[str], cycler.Cycler, str]),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
    )


    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["dpt"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = dpt_args,

    function = dpt_f,

    plot = dpt_plot,

)
