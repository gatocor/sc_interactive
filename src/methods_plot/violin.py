
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

def violin_plot():

    kwargs = config.active_plot_parameters
    fig,ax = plt.subplots()
    
    sc.pl.violin(
        config.adata,
        keys=type_formater(kwargs["keys"],typing.Union[str, typing.Sequence[str]]),
        groupby=type_formater(kwargs["groupby"],typing.Optional[str]),
        log=type_formater(kwargs["log"],bool),
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        stripplot=type_formater(kwargs["stripplot"],bool),
        jitter=type_formater(kwargs["jitter"],typing.Union[float, bool]),
        size=type_formater(kwargs["size"],int),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        scale=type_formater(kwargs["scale"],typing.Literal['area', 'count', 'width']),
        order=type_formater(kwargs["order"],typing.Optional[typing.Sequence[str]]),
        multi_panel=type_formater(kwargs["multi_panel"],typing.Optional[bool]),
        xlabel=type_formater(kwargs["xlabel"],str),
        ylabel=type_formater(kwargs["ylabel"],typing.Union[str, typing.Sequence[str], str]),
        rotation=type_formater(kwargs["rotation"],typing.Optional[float]),
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

config.methods_plot["violin"] = dict(
    
    args = [
    dict(
        input='Input', 
        name='keys', 
        description="typing.Union[str, typing.Sequence[str]]", 
        visible=dict(function="True or config.show_plot"),
        properties=dict(value="''",type="text")
    ),
    dict(
        input='Input', 
        name='groupby', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['groupby'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='log', 
        description="<class 'bool'>", 
        visible=dict(function="str(False)!=config.active_plot_parameters['log'] or config.show_plot"),
        properties=dict(value="False",type="text")
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
        name='stripplot', 
        description="<class 'bool'>", 
        visible=dict(function="str(True)!=config.active_plot_parameters['stripplot'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='jitter', 
        description="typing.Union[float, bool]", 
        visible=dict(function="str(True)!=config.active_plot_parameters['jitter'] or config.show_plot"),
        properties=dict(value="True",type="text")
    ),
    dict(
        input='Input', 
        name='size', 
        description="<class 'int'>", 
        visible=dict(function="str(1)!=config.active_plot_parameters['size'] or config.show_plot"),
        properties=dict(value="1",type="text")
    ),
    dict(
        input='Input', 
        name='layer', 
        description="typing.Optional[str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['layer'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='scale', 
        description="typing.Literal['area', 'count', 'width']", 
        visible=dict(function="'width'!=eval(config.active_plot_parameters['scale']) or config.show_plot"),
        properties=dict(value="'width'",type="text")
    ),
    dict(
        input='Input', 
        name='order', 
        description="typing.Optional[typing.Sequence[str]]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['order'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='multi_panel', 
        description="typing.Optional[bool]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['multi_panel'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='xlabel', 
        description="<class 'str'>", 
        visible=dict(function="''!=eval(config.active_plot_parameters['xlabel']) or config.show_plot"),
        properties=dict(value="''",type="text")
    ),
    dict(
        input='Input', 
        name='ylabel', 
        description="typing.Union[str, typing.Sequence[str], str]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ylabel'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='rotation', 
        description="typing.Optional[float]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['rotation'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),
    dict(
        input='Input', 
        name='ax', 
        description="typing.Optional[matplotlib.axes._axes.Axes]", 
        visible=dict(function="str(None)!=config.active_plot_parameters['ax'] or config.show_plot"),
        properties=dict(value="None",type="text")
    ),],

    function = violin_plot
)
