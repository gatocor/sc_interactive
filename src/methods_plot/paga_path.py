
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

def paga_path_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.paga_path(
        config.adata,
        use_raw=type_formater(kwargs["use_raw"],bool),
        annotations=type_formater(kwargs["annotations"],typing.Sequence[str]),
        color_map=type_formater(kwargs["color_map"],typing.Union[str, matplotlib.colors.Colormap, str]),
        color_maps_annotations=type_formater(kwargs["color_maps_annotations"],typing.Mapping[str, typing.Union[str, matplotlib.colors.Colormap]]),
        palette_groups=type_formater(kwargs["palette_groups"],typing.Optional[typing.Sequence[str]]),
        n_avg=type_formater(kwargs["n_avg"],int),
        groups_key=type_formater(kwargs["groups_key"],typing.Optional[str]),
        xlim=type_formater(kwargs["xlim"],typing.Tuple[typing.Optional[int], typing.Optional[int]]),
        title=type_formater(kwargs["title"],typing.Optional[str]),
        ytick_fontsize=type_formater(kwargs["ytick_fontsize"],typing.Optional[int]),
        title_fontsize=type_formater(kwargs["title_fontsize"],typing.Optional[int]),
        show_node_names=type_formater(kwargs["show_node_names"],bool),
        show_yticks=type_formater(kwargs["show_yticks"],bool),
        show_colorbar=type_formater(kwargs["show_colorbar"],bool),
        legend_fontsize=type_formater(kwargs["legend_fontsize"],typing.Union[int, float, typing.Literal['xx-small', 'x-small', 'small', 'medium', 'large', 'x-large', 'xx-large'], str]),
        legend_fontweight=type_formater(kwargs["legend_fontweight"],typing.Union[int, typing.Literal['light', 'normal', 'medium', 'semibold', 'bold', 'heavy', 'black'], str]),
        normalize_to_zero_one=type_formater(kwargs["normalize_to_zero_one"],bool),
        as_heatmap=type_formater(kwargs["as_heatmap"],bool),
        return_data=type_formater(kwargs["return_data"],bool),
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

config.methods_plot["paga_path"] = dict(
    
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

    function = paga_path_plot
)
