
import numpy
from numpy import inf
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

def scatter_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.scatter(
        config.adata,
        x=type_formater(kwargs["x"],typing.Optional[str]),
        y=type_formater(kwargs["y"],typing.Optional[str]),
        color=type_formater(kwargs["color"],typing.Union[str, typing.Collection[str]]),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        layers=type_formater(kwargs["layers"],typing.Union[str, typing.Collection[str]]),
        sort_order=type_formater(kwargs["sort_order"],bool),
        alpha=type_formater(kwargs["alpha"],typing.Optional[float]),
        basis=type_formater(kwargs["basis"],typing.Optional[typing.Literal['pca', 'tsne', 'umap', 'diffmap', 'draw_graph_fr']]),
        groups=type_formater(kwargs["groups"],typing.Union[str, typing.Iterable[str]]),
        components=type_formater(kwargs["components"],typing.Union[str, typing.Collection[str]]),
        projection=type_formater(kwargs["projection"],typing.Literal['2d', '3d']),
        legend_loc=type_formater(kwargs["legend_loc"],str),
        legend_fontsize=type_formater(kwargs["legend_fontsize"],typing.Union[int, float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], str]),
        legend_fontweight=type_formater(kwargs["legend_fontweight"],typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'], str]),
        legend_fontoutline=type_formater(kwargs["legend_fontoutline"],float),
        color_map=type_formater(kwargs["color_map"],typing.Union[str, matplotlib.colors.Colormap]),
        palette=type_formater(kwargs["palette"],typing.Union[cycler.Cycler, matplotlib.colors.ListedColormap, str, typing.Tuple[float, ...], typing.Sequence[typing.Union[str, typing.Tuple[float, ...]]]]),
        frameon=type_formater(kwargs["frameon"],typing.Optional[bool]),
        right_margin=type_formater(kwargs["right_margin"],typing.Optional[float]),
        left_margin=type_formater(kwargs["left_margin"],typing.Optional[float]),
        size=type_formater(kwargs["size"],typing.Union[int, float, str]),
        title=type_formater(kwargs["title"],typing.Optional[str]),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[str, bool, str]),
        ax=type_formater(kwargs["ax"],typing.Optional[matplotlib.axes._axes.Axes]),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["scatter"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='min_dist', 
        description="<class 'float'>", 
        visible=dict(function="str(0.5)!=config.active_node_parameters['min_dist'] or config.show_parameters"),
        properties=dict(value="0.5",type="text")
    ),
    dict(
        input='Input', 
        name='spread', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['spread'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='n_components', 
        description="<class 'int'>", 
        visible=dict(function="str(2)!=config.active_node_parameters['n_components'] or config.show_parameters"),
        properties=dict(value="2",type="text")
    ),
    dict(
        input='Input', 
        name='maxiter', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_node_parameters['maxiter'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='alpha', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['alpha'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='gamma', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_node_parameters['gamma'] or config.show_parameters"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='negative_sample_rate', 
        description="<class 'int'>", 
        visible=dict(function="str(5)!=config.active_node_parameters['negative_sample_rate'] or config.show_parameters"),
        properties=dict(value="5",type="text")
    ),
    dict(
        input='Input', 
        name='init_pos', 
        description="typing.Union[typing.Literal['paga', 'spectral', 'random'], numpy.ndarray, str]", 
        visible=dict(function="'spectral'!=eval(config.active_node_parameters['init_pos']) or config.show_parameters"),
        properties=dict(value="'spectral'",type="text")
    ),
    dict(
        input='Input', 
        name='random_state', 
        description="typing.Union[str, int, numpy.random.mtrand.RandomState]", 
        visible=dict(function="str(0)!=config.active_node_parameters['random_state'] or config.show_parameters"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='a', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['a'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='b', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_node_parameters['b'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='copy', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_node_parameters['copy'] or config.show_parameters"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='method', 
        description="typing.Literal['umap', 'rapids']", 
        visible=dict(function="'umap'!=eval(config.active_node_parameters['method']) or config.show_parameters"),
        properties=dict(value="'umap'",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_node_parameters['neighbors_key'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),],

    function = scatter_plot
)
