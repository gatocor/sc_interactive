
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

neighbors_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='n_neighbors', 
        description="<class 'int'>", 
        visible=dict(function="str(15)!=get_node(config.selected)['data']['parameters']['n_neighbors'] or config.show_parameters"),
        properties=dict(value="15",type="text")
    ),
    dict(
        input='Input', 
        name='n_pcs', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['n_pcs'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_rep', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['use_rep'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='knn', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['parameters']['knn'] or config.show_parameters"),
        properties=dict(value="True",type="text")
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
        name='method', 
        description="typing.Optional[typing.Literal['umap', 'gauss', 'rapids']]", 
        visible=dict(function="'umap'!=eval(get_node(config.selected)['data']['parameters']['method']) or config.show_parameters"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='metric', 
        description="typing.Union[typing.Literal['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'], typing.Literal['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'], typing.Callable[[numpy.ndarray, numpy.ndarray], float]]", 
        visible=dict(function="'euclidean'!=eval(get_node(config.selected)['data']['parameters']['metric']) or config.show_parameters"),
        properties=dict(value="'euclidean'",type="text")
    ),
    dict(
        input='Input', 
        name='metric_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=get_node(config.selected)['data']['parameters']['metric_kwds'] or config.show_parameters"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='key_added', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['key_added'] or config.show_parameters"),
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

def neighbors_f(adata,kwargs):

    sc.pp.neighbors(
        adata,
        n_neighbors=type_formater(kwargs["n_neighbors"],int),
        n_pcs=type_formater(kwargs["n_pcs"],typing.Optional[int]),
        use_rep=type_formater(kwargs["use_rep"],typing.Optional[str]),
        knn=type_formater(kwargs["knn"],bool),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        method=type_formater(kwargs["method"],typing.Optional[typing.Literal['umap', 'gauss', 'rapids']]),
        metric=type_formater(kwargs["metric"],typing.Union[typing.Literal['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan'], typing.Literal['braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'], typing.Callable[[numpy.ndarray, numpy.ndarray], float]]),
        metric_kwds=type_formater(kwargs["metric_kwds"],typing.Mapping[str, typing.Any]),
        key_added=type_formater(kwargs["key_added"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def neighbors_plot():

    kwargs = get_node(config.selected)['data']['plot']
    

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["neighbors"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = neighbors_args,

    function = neighbors_f,

    plot = neighbors_plot,

)
