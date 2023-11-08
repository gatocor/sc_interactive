
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

def paga_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.paga(
        config.adata,
        threshold=type_formater(kwargs["threshold"],typing.Optional[float]),
        color=type_formater(kwargs["color"],typing.Union[str, typing.Mapping[typing.Union[str, int], typing.Mapping[typing.Any, float]], str]),
        layout=type_formater(kwargs["layout"],typing.Optional[typing.Literal['fa', 'fr', 'rt', 'rt_circular', 'drl', 'eq_tree', ...]]),
        layout_kwds=type_formater(kwargs["layout_kwds"],typing.Mapping[str, typing.Any]),
        init_pos=type_formater(kwargs["init_pos"],typing.Optional[numpy.ndarray]),
        root=type_formater(kwargs["root"],typing.Union[int, str, typing.Sequence[int], str]),
        labels=type_formater(kwargs["labels"],typing.Union[str, typing.Sequence[str], typing.Mapping[str, str], str]),
        single_component=type_formater(kwargs["single_component"],bool),
        solid_edges=type_formater(kwargs["solid_edges"],str),
        dashed_edges=type_formater(kwargs["dashed_edges"],typing.Optional[str]),
        transitions=type_formater(kwargs["transitions"],typing.Optional[str]),
        fontsize=type_formater(kwargs["fontsize"],typing.Optional[int]),
        fontweight=type_formater(kwargs["fontweight"],str),
        fontoutline=type_formater(kwargs["fontoutline"],typing.Optional[int]),
        text_kwds=type_formater(kwargs["text_kwds"],typing.Mapping[str, typing.Any]),
        node_size_scale=type_formater(kwargs["node_size_scale"],float),
        node_size_power=type_formater(kwargs["node_size_power"],float),
        edge_width_scale=type_formater(kwargs["edge_width_scale"],float),
        min_edge_width=type_formater(kwargs["min_edge_width"],typing.Optional[float]),
        max_edge_width=type_formater(kwargs["max_edge_width"],typing.Optional[float]),
        arrowsize=type_formater(kwargs["arrowsize"],int),
        title=type_formater(kwargs["title"],typing.Optional[str]),
        left_margin=type_formater(kwargs["left_margin"],float),
        random_state=type_formater(kwargs["random_state"],typing.Optional[int]),
        pos=type_formater(kwargs["pos"],typing.Union[numpy.ndarray, str, pathlib.Path, str]),
        normalize_to_color=type_formater(kwargs["normalize_to_color"],bool),
        cmap=type_formater(kwargs["cmap"],typing.Union[str, matplotlib.colors.Colormap]),
        cax=type_formater(kwargs["cax"],typing.Optional[matplotlib.axes._axes.Axes]),
        cb_kwds=type_formater(kwargs["cb_kwds"],typing.Mapping[str, typing.Any]),
        frameon=type_formater(kwargs["frameon"],typing.Optional[bool]),
        add_pos=type_formater(kwargs["add_pos"],bool),
        export_to_gexf=type_formater(kwargs["export_to_gexf"],bool),
        use_raw=type_formater(kwargs["use_raw"],bool),
        plot=type_formater(kwargs["plot"],bool),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
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

config.methods_plot["paga"] = dict(
    
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

    function = paga_plot
)
