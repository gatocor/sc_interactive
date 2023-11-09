
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

def paga_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.paga(
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
        cax=ax,
        cb_kwds=type_formater(kwargs["cb_kwds"],typing.Mapping[str, typing.Any]),
        frameon=type_formater(kwargs["frameon"],typing.Optional[bool]),
        add_pos=type_formater(kwargs["add_pos"],bool),
        export_to_gexf=type_formater(kwargs["export_to_gexf"],bool),
        use_raw=type_formater(kwargs["use_raw"],bool),
        plot=type_formater(kwargs["plot"],bool),
        ax=ax,
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png", transparent=True)
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["paga"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='threshold', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['threshold'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='color', 
        description="typing.Union[str, typing.Mapping[typing.Union[str, int], typing.Mapping[typing.Any, float]], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layout', 
        description="typing.Optional[typing.Literal['fa', 'fr', 'rt', 'rt_circular', 'drl', 'eq_tree', ...]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['layout'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layout_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=config.active_plot_parameters['layout_kwds'] or config.show_plot"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='init_pos', 
        description="typing.Optional[numpy.ndarray]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['init_pos'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='root', 
        description="typing.Union[int, str, typing.Sequence[int], str]", 
        visible=dict(function="str(0)!=config.active_plot_parameters['root'] or config.show_plot"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='labels', 
        description="typing.Union[str, typing.Sequence[str], typing.Mapping[str, str], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['labels'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='single_component', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['single_component'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='solid_edges', 
        description="<class 'str'>", 
        visible=dict(function="'connectivities'!=eval(config.active_plot_parameters['solid_edges']) or config.show_plot"),
        properties=dict(value="'connectivities'",type="text")
    ),
    dict(
        input='Input', 
        name='dashed_edges', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['dashed_edges'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='transitions', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['transitions'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='fontsize', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['fontsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='fontweight', 
        description="<class 'str'>", 
        visible=dict(function="'bold'!=eval(config.active_plot_parameters['fontweight']) or config.show_plot"),
        properties=dict(value="'bold'",type="text")
    ),
    dict(
        input='Input', 
        name='fontoutline', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['fontoutline'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='text_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=config.active_plot_parameters['text_kwds'] or config.show_plot"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='node_size_scale', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_plot_parameters['node_size_scale'] or config.show_plot"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='node_size_power', 
        description="<class 'float'>", 
        visible=dict(function="str(0.5)!=config.active_plot_parameters['node_size_power'] or config.show_plot"),
        properties=dict(value="0.5",type="text")
    ),
    dict(
        input='Input', 
        name='edge_width_scale', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=config.active_plot_parameters['edge_width_scale'] or config.show_plot"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='min_edge_width', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['min_edge_width'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_edge_width', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['max_edge_width'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='arrowsize', 
        description="<class 'int'>", 
        visible=dict(function="str(30)!=config.active_plot_parameters['arrowsize'] or config.show_plot"),
        properties=dict(value="30",type="text")
    ),
    dict(
        input='Input', 
        name='title', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['title'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='left_margin', 
        description="<class 'float'>", 
        visible=dict(function="str(0.01)!=config.active_plot_parameters['left_margin'] or config.show_plot"),
        properties=dict(value="0.01",type="text")
    ),
    dict(
        input='Input', 
        name='random_state', 
        description="typing.Optional[int]", 
        visible=dict(function="str(0)!=config.active_plot_parameters['random_state'] or config.show_plot"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='pos', 
        description="typing.Union[numpy.ndarray, str, pathlib.Path, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['pos'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='normalize_to_color', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['normalize_to_color'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='cmap', 
        description="typing.Union[str, matplotlib.colors.Colormap]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['cmap'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['cax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cb_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=config.active_plot_parameters['cb_kwds'] or config.show_plot"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='frameon', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['frameon'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='add_pos', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['add_pos'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='export_to_gexf', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['export_to_gexf'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='use_raw', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['use_raw'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='plot', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['plot'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = paga_plot
)
