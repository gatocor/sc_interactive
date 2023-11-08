
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

def scatter_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.scatter(
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

config.methods_plot["scatter"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='x', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['x'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='y', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['y'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='color', 
        description="typing.Union[str, typing.Collection[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color'] or config.show_plot"),
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
        name='layers', 
        description="typing.Union[str, typing.Collection[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['layers'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='sort_order', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['sort_order'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='alpha', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['alpha'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='basis', 
        description="typing.Optional[typing.Literal['pca', 'tsne', 'umap', 'diffmap', 'draw_graph_fr']]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['basis'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='groups', 
        description="typing.Union[str, typing.Iterable[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groups'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='components', 
        description="typing.Union[str, typing.Collection[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['components'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='projection', 
        description="typing.Literal['2d', '3d']", 
        visible=dict(function="'2d'!=eval(config.active_plot_parameters['projection']) or config.show_plot"),
        properties=dict(value="'2d'",type="text")
    ),
    dict(
        input='Input', 
        name='legend_loc', 
        description="<class 'str'>", 
        visible=dict(function="'right margin'!=eval(config.active_plot_parameters['legend_loc']) or config.show_plot"),
        properties=dict(value="'right margin'",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontsize', 
        description="typing.Union[int, float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['legend_fontsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontweight', 
        description="typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['legend_fontweight'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='legend_fontoutline', 
        description="<class 'float'>", 
        visible=dict(function="str(None)!=config.active_plot_parameters['legend_fontoutline'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[str, matplotlib.colors.Colormap]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color_map'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='palette', 
        description="typing.Union[cycler.Cycler, matplotlib.colors.ListedColormap, str, typing.Tuple[float, ...], typing.Sequence[typing.Union[str, typing.Tuple[float, ...]]]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['palette'] or config.show_plot"),
        properties=dict(value="None",type="text")
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
        name='right_margin', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['right_margin'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='left_margin', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['left_margin'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='size', 
        description="typing.Union[int, float, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['size'] or config.show_plot"),
        properties=dict(value="None",type="text")
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
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = scatter_plot
)
