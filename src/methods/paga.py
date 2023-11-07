
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

paga_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='groups', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['groups']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_rna_velocity', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['parameters']['use_rna_velocity']"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='model', 
        description="typing.Literal['v1.2', 'v1.0']", 
        visible=dict(function="'v1.2'!=eval(get_node(config.selected)['data']['parameters']['model'])"),
        properties=dict(value="'v1.2'",type="text")
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
        name='threshold', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['threshold']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='color', 
        description="typing.Union[str, typing.Mapping[typing.Union[str, int], typing.Mapping[typing.Any, float]], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['color']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layout', 
        description="typing.Optional[typing.Literal['fa', 'fr', 'rt', 'rt_circular', 'drl', 'eq_tree', ...]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['layout']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layout_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=get_node(config.selected)['data']['plot']['layout_kwds']"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='init_pos', 
        description="typing.Optional[numpy.ndarray]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['init_pos']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='root', 
        description="typing.Union[int, str, typing.Sequence[int], str]", 
        visible=dict(function="str(0)!=get_node(config.selected)['data']['plot']['root']"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='labels', 
        description="typing.Union[str, typing.Sequence[str], typing.Mapping[str, str], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['labels']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='single_component', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['plot']['single_component']"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='solid_edges', 
        description="<class 'str'>", 
        visible=dict(function="'connectivities'!=eval(get_node(config.selected)['data']['plot']['solid_edges'])"),
        properties=dict(value="'connectivities'",type="text")
    ),
    dict(
        input='Input', 
        name='dashed_edges', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['dashed_edges']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='transitions', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['transitions']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='fontsize', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['fontsize']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='fontweight', 
        description="<class 'str'>", 
        visible=dict(function="'bold'!=eval(get_node(config.selected)['data']['plot']['fontweight'])"),
        properties=dict(value="'bold'",type="text")
    ),
    dict(
        input='Input', 
        name='fontoutline', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['fontoutline']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='text_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=get_node(config.selected)['data']['plot']['text_kwds']"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='node_size_scale', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=get_node(config.selected)['data']['plot']['node_size_scale']"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='node_size_power', 
        description="<class 'float'>", 
        visible=dict(function="str(0.5)!=get_node(config.selected)['data']['plot']['node_size_power']"),
        properties=dict(value="0.5",type="text")
    ),
    dict(
        input='Input', 
        name='edge_width_scale', 
        description="<class 'float'>", 
        visible=dict(function="str(1.0)!=get_node(config.selected)['data']['plot']['edge_width_scale']"),
        properties=dict(value="1.0",type="text")
    ),
    dict(
        input='Input', 
        name='min_edge_width', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['min_edge_width']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='max_edge_width', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['max_edge_width']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='arrowsize', 
        description="<class 'int'>", 
        visible=dict(function="str(30)!=get_node(config.selected)['data']['plot']['arrowsize']"),
        properties=dict(value="30",type="text")
    ),
    dict(
        input='Input', 
        name='title', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['title']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='left_margin', 
        description="<class 'float'>", 
        visible=dict(function="str(0.01)!=get_node(config.selected)['data']['plot']['left_margin']"),
        properties=dict(value="0.01",type="text")
    ),
    dict(
        input='Input', 
        name='random_state', 
        description="typing.Optional[int]", 
        visible=dict(function="str(0)!=get_node(config.selected)['data']['plot']['random_state']"),
        properties=dict(value="0",type="text")
    ),
    dict(
        input='Input', 
        name='pos', 
        description="typing.Union[numpy.ndarray, str, pathlib.Path, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['pos']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='normalize_to_color', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['plot']['normalize_to_color']"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='cmap', 
        description="typing.Union[str, matplotlib.colors.Colormap]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['cmap']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['cax']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cb_kwds', 
        description="typing.Mapping[str, typing.Any]", 
        visible=dict(function="str({})!=get_node(config.selected)['data']['plot']['cb_kwds']"),
        properties=dict(value="{}",type="text")
    ),
    dict(
        input='Input', 
        name='frameon', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['frameon']"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='add_pos', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['plot']['add_pos']"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='export_to_gexf', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['plot']['export_to_gexf']"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='use_raw', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['plot']['use_raw']"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='plot', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['plot']['plot']"),
        properties=dict(value="True",type="text")
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
    ),]
)

def paga_f(adata,kwargs):

    sc.tl.paga(
        adata,
        groups=type_formater(kwargs["groups"],typing.Optional[str]),
        use_rna_velocity=type_formater(kwargs["use_rna_velocity"],bool),
        model=type_formater(kwargs["model"],typing.Literal['v1.2', 'v1.0']),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def paga_plot():

    kwargs = get_node(config.selected)['data']['plot']
    
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

config.methods["paga"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = paga_args,

    function = paga_f,

    plot = paga_plot,

)
