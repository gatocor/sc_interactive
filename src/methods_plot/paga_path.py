
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

def paga_path_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.paga_path(
        config.adata,
        nodes=type_formater(kwargs["nodes"],typing.Sequence[typing.Union[str, int]]),
        keys=type_formater(kwargs["keys"],typing.Sequence[str]),
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

config.methods_plot["paga_path"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='nodes', 
        description="typing.Sequence[typing.Union[str, int]]", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
    ),
    dict(
        input='Input', 
        name='keys', 
        description="typing.Sequence[str]", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
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
        name='annotations', 
        description="typing.Sequence[str]", 
        visible=dict(function="str(('dpt_pseudotime',))!=config.active_plot_parameters['annotations'] or config.show_plot"),
        properties=dict(value="('dpt_pseudotime',)",type="text")
    ),
    dict(
        input='Input', 
        name='color_map', 
        description="typing.Union[str, matplotlib.colors.Colormap, str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['color_map'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='color_maps_annotations', 
        description="typing.Mapping[str, typing.Union[str, matplotlib.colors.Colormap]]", 
        visible=dict(function="str({'dpt_pseudotime': 'Greys'})!=config.active_plot_parameters['color_maps_annotations'] or config.show_plot"),
        properties=dict(value="{'dpt_pseudotime': 'Greys'}",type="text")
    ),
    dict(
        input='Input', 
        name='palette_groups', 
        description="typing.Optional[typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['palette_groups'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='n_avg', 
        description="<class 'int'>", 
        visible=dict(function="str(1)!=config.active_plot_parameters['n_avg'] or config.show_plot"),
        properties=dict(value="1",type="text")
    ),
    dict(
        input='Input', 
        name='groups_key', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groups_key'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='xlim', 
        description="typing.Tuple[typing.Optional[int], typing.Optional[int]]", 
        visible=dict(function="str((None, None))!=config.active_plot_parameters['xlim'] or config.show_plot"),
        properties=dict(value="(None, None)",type="text")
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
        name='ytick_fontsize', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ytick_fontsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='title_fontsize', 
        description="typing.Optional[int]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['title_fontsize'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='show_node_names', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['show_node_names'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='show_yticks', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['show_yticks'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='show_colorbar', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['show_colorbar'] or config.show_plot"),
        properties=dict(value="True",type="text")
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
        name='normalize_to_zero_one', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['normalize_to_zero_one'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='as_heatmap', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['as_heatmap'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='return_data', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['return_data'] or config.show_plot"),
        properties=dict(value="False",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = paga_path_plot
)
