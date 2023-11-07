
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

louvain_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='resolution', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['resolution'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='random_state', 
        description="typing.Union[str, int, numpy.random.mtrand.RandomState]", 
        visible=dict(function="str(0)!=get_node(config.selected)['data']['parameters']['random_state'] or config.show_parameters"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='restrict_to', 
        description="typing.Optional[typing.Tuple[str, typing.Sequence[str]]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['restrict_to'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='key_added', 
        description="<class 'str'>", 
        visible=dict(function="'louvain'!=eval(get_node(config.selected)['data']['parameters']['key_added']) or config.show_parameters"),
        properties=dict(value="'louvain'",type="text")
    ),
    dict(
        input='Input', 
        name='adjacency', 
        description="typing.Optional[scipy.sparse._matrix.spmatrix]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['adjacency'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='flavor', 
        description="typing.Literal['vtraag', 'igraph', 'rapids']", 
        visible=dict(function="'vtraag'!=eval(get_node(config.selected)['data']['parameters']['flavor']) or config.show_parameters"),
        properties=dict(value="'vtraag'",type="text")
    ),
    dict(
        input='Input', 
        name='directed', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['directed'] or config.show_parameters"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='use_weights', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['use_weights'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='partition_type', 
        description="typing.Optional[typing.Type[louvain.VertexPartition.MutableVertexPartition]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['partition_type'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='partition_kwargs', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=get_node(config.selected)['data']['parameters']['partition_kwargs'] or config.show_parameters"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['neighbors_key'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='obsp', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['obsp'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),],
    postexecution = [],
    plot = []
)

def louvain_f(adata,kwargs):

    sc.tl.louvain(
        adata,
        resolution=type_formater(kwargs["resolution"],typing.Optional[float]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        restrict_to=type_formater(kwargs["restrict_to"],typing.Optional[typing.Tuple[str, typing.Sequence[str]]]),
        key_added=type_formater(kwargs["key_added"],str),
        adjacency=type_formater(kwargs["adjacency"],typing.Optional[scipy.sparse._matrix.spmatrix]),
        flavor=type_formater(kwargs["flavor"],typing.Literal['vtraag', 'igraph', 'rapids']),
        directed=type_formater(kwargs["directed"],bool),
        use_weights=type_formater(kwargs["use_weights"],bool),
        partition_type=type_formater(kwargs["partition_type"],typing.Optional[typing.Type[louvain.VertexPartition.MutableVertexPartition]]),
        partition_kwargs=type_formater(kwargs["partition_kwargs"],typing.Mapping[str, typing.Any]),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        obsp=type_formater(kwargs["obsp"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def louvain_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["louvain"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = louvain_args,

    function = louvain_f,

    plot = louvain_plot,

)
