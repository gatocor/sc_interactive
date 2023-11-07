
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

embedding_density_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='basis', 
        description="<class 'str'>", 
        visible=dict(function="'umap'!=eval(get_node(config.selected)['data']['parameters']['basis'])"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='groupby', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['groupby']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='key_added', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['key_added']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='components', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['components']"),
        properties=dict(value="None",type="text")
    ),],
    postexecution = [],
    plot = [
    dict(
        input='Input', 
        name='basis', 
        description="<class 'str'>", 
        visible=dict(function="'umap'!=eval(get_node(config.selected)['data']['plot']['basis'])"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['key']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='groupby', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['groupby']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='group', 
        description="typing.Union[str, typing.List[str], str]", 
        visible=dict(function="'all'!=eval(get_node(config.selected)['data']['plot']['group'])"),
        properties=dict(value="'all'",type="text")
    ),
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[matplotlib.colors.Colormap, str]", 
        visible=dict(function="'YlOrRd'!=eval(get_node(config.selected)['data']['plot']['color_map'])"),
        properties=dict(value="'YlOrRd'",type="text")
    ),
    dict(
        input='Input', 
        name='bg_dotsize', 
        description="typing.Optional[int]", 
        visible=dict(function="str(80)!=get_node(config.selected)['data']['plot']['bg_dotsize']"),
        properties=dict(value="80",type="text")
    ),
    dict(
        input='Input', 
        name='fg_dotsize', 
        description="typing.Optional[int]", 
        visible=dict(function="str(180)!=get_node(config.selected)['data']['plot']['fg_dotsize']"),
        properties=dict(value="180",type="text")
    ),
    dict(
        input='Input', 
        name='vmax', 
        description="typing.Optional[int]", 
        visible=dict(function="str(1)!=get_node(config.selected)['data']['plot']['vmax']"),
        properties=dict(value="1",type="text")
    ),
    dict(
        input='Input', 
        name='vmin', 
        description="typing.Optional[int]", 
        visible=dict(function="str(0)!=get_node(config.selected)['data']['plot']['vmin']"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='vcenter', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['vcenter']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='norm', 
        description="typing.Optional[matplotlib.colors.Normalize]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['norm']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='ncols', 
        description="typing.Optional[int]", 
        visible=dict(function="str(4)!=get_node(config.selected)['data']['plot']['ncols']"),
        properties=dict(value="4",type="text")
    ),
    dict(
        input='Input', 
        name='hspace', 
        description="typing.Optional[float]", 
        visible=dict(function="str(0.25)!=get_node(config.selected)['data']['plot']['hspace']"),
        properties=dict(value="0.25",type="text")
    ),
    dict(
        input='Input', 
        name='wspace', 
        description="<class 'str'>", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['wspace']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='title', 
        description="<class 'str'>", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['title']"),
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
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['ax']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='return_fig', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['return_fig']"),
        properties=dict(value="None",type="text")
    ),]
)

def embedding_density_f(adata,kwargs):

    sc.tl.embedding_density(
        adata,
        basis=type_formater(kwargs["basis"],str),
        groupby=type_formater(kwargs["groupby"],typing.Optional[str]),
        key_added=type_formater(kwargs["key_added"],typing.Optional[str]),
        components=type_formater(kwargs["components"],typing.Union[str, typing.Sequence[str]]),
    )
        
    return

def embedding_density_plot():

    kwargs = get_node(config.selected)['data']['plot']
    
    fig = sc.pl.embedding_density(
        config.adata,
        basis=type_formater(kwargs["basis"],str),
        key=type_formater(kwargs["key"],typing.Optional[str]),
        groupby=type_formater(kwargs["groupby"],typing.Optional[str]),
        group=type_formater(kwargs["group"],typing.Union[str, typing.List[str], str]),
        color_map=type_formater(kwargs["color_map"],typing.Union[matplotlib.colors.Colormap, str]),
        bg_dotsize=type_formater(kwargs["bg_dotsize"],typing.Optional[int]),
        fg_dotsize=type_formater(kwargs["fg_dotsize"],typing.Optional[int]),
        vmax=type_formater(kwargs["vmax"],typing.Optional[int]),
        vmin=type_formater(kwargs["vmin"],typing.Optional[int]),
        vcenter=type_formater(kwargs["vcenter"],typing.Optional[int]),
        norm=type_formater(kwargs["norm"],typing.Optional[matplotlib.colors.Normalize]),
        ncols=type_formater(kwargs["ncols"],typing.Optional[int]),
        hspace=type_formater(kwargs["hspace"],typing.Optional[float]),
        wspace=type_formater(kwargs["wspace"],str),
        title=type_formater(kwargs["title"],str),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
        ax=type_formater(kwargs["ax"],typing.Optional[matplotlib.axes._axes.Axes]),
        return_fig=True,
    )


    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["embedding_density"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = embedding_density_args,

    function = embedding_density_f,

    plot = embedding_density_plot,

)
