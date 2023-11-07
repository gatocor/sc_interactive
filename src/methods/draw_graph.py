
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

draw_graph_args = dict(
    execution = [ARGINPUT,
    dict(
        input='Input', 
        name='layout', 
        description="typing.Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa']", 
        visible=dict(function="'fa'!=eval(get_node(config.selected)['data']['parameters']['layout']) or config.show_parameters"),
        properties=dict(value="'fa'",type="text")
    ),
    dict(
        input='Input', 
        name='init_pos', 
        description="typing.Union[str, bool, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['init_pos'] or config.show_parameters"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='root', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['root'] or config.show_parameters"),
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
        name='n_jobs', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['n_jobs'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
        name='key_added_ext', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['parameters']['key_added_ext'] or config.show_parameters"),
        properties=dict(value="None",type="text")
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
    plot = [
    dict(
        input='Input', 
        name='color', 
        description="typing.Union[str, typing.Sequence[str], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['color'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='gene_symbols', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['gene_symbols'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='use_raw', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['use_raw'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='sort_order', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['plot']['sort_order'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='edges', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['plot']['edges'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='edges_width', 
        description="<class 'float'>", 
        visible=dict(function="str(0.1)!=get_node(config.selected)['data']['plot']['edges_width'] or config.show_plot"),
        properties=dict(value="0.1",type="text")
    ),
    dict(
        input='Input', 
        name='edges_color', 
        description="typing.Union[str, typing.Sequence[float], typing.Sequence[str]]", 
        visible=dict(function="'grey'!=eval(get_node(config.selected)['data']['plot']['edges_color']) or config.show_plot"),
        properties=dict(value="'grey'",type="text")
    ),
    dict(
        input='Input', 
        name='neighbors_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['neighbors_key'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='arrows', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['plot']['arrows'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='arrows_kwds', 
        description="typing.Optional[typing.Mapping[str, typing.Any]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['arrows_kwds'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='groups', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['groups'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='components', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['components'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='dimensions', 
        description="typing.Union[typing.Tuple[int, int], typing.Sequence[typing.Tuple[int, int]], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['dimensions'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['layer'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='projection', 
        description="typing.Literal['2d', '3d']", 
        visible=dict(function="'2d'!=eval(get_node(config.selected)['data']['plot']['projection']) or config.show_plot"),
        properties=dict(value="'2d'",type="text")
    ),
    dict(
        input='Input', 
        name='scale_factor', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['scale_factor'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[matplotlib.colors.Colormap, str, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['color_map'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='cmap', 
        description="typing.Union[matplotlib.colors.Colormap, str, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['cmap'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='palette', 
        description="typing.Union[str, typing.Sequence[str], cycler.Cycler, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['palette'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='na_color', 
        description="typing.Union[str, typing.Tuple[float, ...]]", 
        visible=dict(function="'lightgray'!=eval(get_node(config.selected)['data']['plot']['na_color']) or config.show_plot"),
        properties=dict(value="'lightgray'",type="text")
    ),
    dict(
        input='Input', 
        name='na_in_legend', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=get_node(config.selected)['data']['plot']['na_in_legend'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='size', 
        description="typing.Union[float, typing.Sequence[float], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['size'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='frameon', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['frameon'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontsize', 
        description="typing.Union[int, float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['legend_fontsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontweight', 
        description="typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]", 
        visible=dict(function="'bold'!=eval(get_node(config.selected)['data']['plot']['legend_fontweight']) or config.show_plot"),
        properties=dict(value="'bold'",type="text")
    ),
    dict(
        input='Input', 
        name='legend_loc', 
        description="<class 'str'>", 
        visible=dict(function="'right margin'!=eval(get_node(config.selected)['data']['plot']['legend_loc']) or config.show_plot"),
        properties=dict(value="'right margin'",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontoutline', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['legend_fontoutline'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='colorbar_loc', 
        description="typing.Optional[str]", 
        visible=dict(function="'right'!=eval(get_node(config.selected)['data']['plot']['colorbar_loc']) or config.show_plot"),
        properties=dict(value="'right'",type="text")
    ),
    dict(
        input='Input', 
        name='vmax', 
        description="typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['vmax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vmin', 
        description="typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['vmin'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='vcenter', 
        description="typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['vcenter'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='norm', 
        description="typing.Union[matplotlib.colors.Normalize, typing.Sequence[matplotlib.colors.Normalize], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['norm'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='add_outline', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(False)!=get_node(config.selected)['data']['plot']['add_outline'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='outline_width', 
        description="typing.Tuple[float, float]", 
        visible=dict(function="str((0.3, 0.05))!=get_node(config.selected)['data']['plot']['outline_width'] or config.show_plot"),
        properties=dict(value="(0.3, 0.05)",type="text")
    ),
    dict(
        input='Input', 
        name='outline_color', 
        description="typing.Tuple[str, str]", 
        visible=dict(function="str(('black', 'white'))!=get_node(config.selected)['data']['plot']['outline_color'] or config.show_plot"),
        properties=dict(value="('black', 'white')",type="text")
    ),
    dict(
        input='Input', 
        name='ncols', 
        description="<class 'int'>", 
        visible=dict(function="str(4)!=get_node(config.selected)['data']['plot']['ncols'] or config.show_plot"),
        properties=dict(value="4",type="text")
    ),
    dict(
        input='Input', 
        name='hspace', 
        description="<class 'float'>", 
        visible=dict(function="str(0.25)!=get_node(config.selected)['data']['plot']['hspace'] or config.show_plot"),
        properties=dict(value="0.25",type="text")
    ),
    dict(
        input='Input', 
        name='wspace', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['wspace'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='title', 
        description="typing.Union[str, typing.Sequence[str], str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['title'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='show', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['show'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='save', 
        description="typing.Union[bool, str, str]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['save'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='return_fig', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['return_fig'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='layout', 
        description="typing.Optional[typing.Literal['fa', 'fr', 'rt', 'rt_circular', 'drl', 'eq_tree', ...]]", 
        visible=dict(function="str(None)!=get_node(config.selected)['data']['plot']['layout'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),]
)

def draw_graph_f(adata,kwargs):

    sc.tl.draw_graph(
        adata,
        layout=type_formater(kwargs["layout"],typing.Literal['fr', 'drl', 'kk', 'grid_fr', 'lgl', 'rt', 'rt_circular', 'fa']),
        init_pos=type_formater(kwargs["init_pos"],typing.Union[str, bool, str]),
        root=type_formater(kwargs["root"],typing.Optional[int]),
        random_state=type_formater(kwargs["random_state"],typing.Union[str, int, numpy.random.mtrand.RandomState]),
        n_jobs=type_formater(kwargs["n_jobs"],typing.Optional[int]),
        adjacency=type_formater(kwargs["adjacency"],typing.Optional[scipy.sparse._matrix.spmatrix]),
        key_added_ext=type_formater(kwargs["key_added_ext"],typing.Optional[str]),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        obsp=type_formater(kwargs["obsp"],typing.Optional[str]),
        copy=type_formater(kwargs["copy"],bool),
    )
        
    return

def draw_graph_plot():

    kwargs = get_node(config.selected)['data']['plot']
    
    fig = sc.pl.draw_graph(
        config.adata,
        color=type_formater(kwargs["color"],typing.Union[str, typing.Sequence[str], str]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        sort_order=type_formater(kwargs["sort_order"],bool),
        edges=type_formater(kwargs["edges"],bool),
        edges_width=type_formater(kwargs["edges_width"],float),
        edges_color=type_formater(kwargs["edges_color"],typing.Union[str, typing.Sequence[float], typing.Sequence[str]]),
        neighbors_key=type_formater(kwargs["neighbors_key"],typing.Optional[str]),
        arrows=type_formater(kwargs["arrows"],bool),
        arrows_kwds=type_formater(kwargs["arrows_kwds"],typing.Optional[typing.Mapping[str, typing.Any]]),
        groups=type_formater(kwargs["groups"],typing.Optional[str]),
        components=type_formater(kwargs["components"],typing.Union[str, typing.Sequence[str]]),
        dimensions=type_formater(kwargs["dimensions"],typing.Union[typing.Tuple[int, int], typing.Sequence[typing.Tuple[int, int]], str]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        projection=type_formater(kwargs["projection"],typing.Literal['2d', '3d']),
        scale_factor=type_formater(kwargs["scale_factor"],typing.Optional[float]),
        color_map=type_formater(kwargs["color_map"],typing.Union[matplotlib.colors.Colormap, str, str]),
        cmap=type_formater(kwargs["cmap"],typing.Union[matplotlib.colors.Colormap, str, str]),
        palette=type_formater(kwargs["palette"],typing.Union[str, typing.Sequence[str], cycler.Cycler, str]),
        na_color=type_formater(kwargs["na_color"],typing.Union[str, typing.Tuple[float, ...]]),
        na_in_legend=type_formater(kwargs["na_in_legend"],bool),
        size=type_formater(kwargs["size"],typing.Union[float, typing.Sequence[float], str]),
        frameon=type_formater(kwargs["frameon"],typing.Optional[bool]),
        legend_fontsize=type_formater(kwargs["legend_fontsize"],typing.Union[int, float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], str]),
        legend_fontweight=type_formater(kwargs["legend_fontweight"],typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black']]),
        legend_loc=type_formater(kwargs["legend_loc"],str),
        legend_fontoutline=type_formater(kwargs["legend_fontoutline"],typing.Optional[int]),
        colorbar_loc=type_formater(kwargs["colorbar_loc"],typing.Optional[str]),
        vmax=type_formater(kwargs["vmax"],typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]),
        vmin=type_formater(kwargs["vmin"],typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]),
        vcenter=type_formater(kwargs["vcenter"],typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float], typing.Sequence[typing.Union[str, float, typing.Callable[[typing.Sequence[float]], float]]], str]),
        norm=type_formater(kwargs["norm"],typing.Union[matplotlib.colors.Normalize, typing.Sequence[matplotlib.colors.Normalize], str]),
        add_outline=type_formater(kwargs["add_outline"],typing.Optional[bool]),
        outline_width=type_formater(kwargs["outline_width"],typing.Tuple[float, float]),
        outline_color=type_formater(kwargs["outline_color"],typing.Tuple[str, str]),
        ncols=type_formater(kwargs["ncols"],int),
        hspace=type_formater(kwargs["hspace"],float),
        wspace=type_formater(kwargs["wspace"],typing.Optional[float]),
        title=type_formater(kwargs["title"],typing.Union[str, typing.Sequence[str], str]),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[bool, str, str]),
        ax=type_formater(kwargs["ax"],typing.Optional[matplotlib.axes._axes.Axes]),
        return_fig=True,
        layout=type_formater(kwargs["layout"],typing.Optional[typing.Literal['fa', 'fr', 'rt', 'rt_circular', 'drl', 'eq_tree', ...]]),
    )


    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods["draw_graph"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = draw_graph_args,

    function = draw_graph_f,

    plot = draw_graph_plot,

)
