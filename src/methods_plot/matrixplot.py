
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

def matrixplot_plot():

    kwargs = config.selected_plot_parameters
    
    fig = sc.pl.matrixplot(
        config.adata,
        use_raw=type_formater(kwargs["use_raw"],typing.Optional[bool]),
        log=type_formater(kwargs["log"],bool),
        num_categories=type_formater(kwargs["num_categories"],int),
        figsize=type_formater(kwargs["figsize"],typing.Optional[typing.Tuple[float, float]]),
        dendrogram=type_formater(kwargs["dendrogram"],typing.Union[bool, str]),
        title=type_formater(kwargs["title"],typing.Optional[str]),
        cmap=type_formater(kwargs["cmap"],typing.Optional[str]),
        colorbar_title=type_formater(kwargs["colorbar_title"],typing.Optional[str]),
        gene_symbols=type_formater(kwargs["gene_symbols"],typing.Optional[str]),
        var_group_positions=type_formater(kwargs["var_group_positions"],typing.Optional[typing.Sequence[typing.Tuple[int, int]]]),
        var_group_labels=type_formater(kwargs["var_group_labels"],typing.Optional[typing.Sequence[str]]),
        var_group_rotation=type_formater(kwargs["var_group_rotation"],typing.Optional[float]),
        layer=type_formater(kwargs["layer"],typing.Optional[str]),
        standard_scale=type_formater(kwargs["standard_scale"],typing.Literal['var', 'group']),
        values_df=type_formater(kwargs["values_df"],typing.Optional[pandas.core.frame.DataFrame]),
        swap_axes=type_formater(kwargs["swap_axes"],bool),
        show=type_formater(kwargs["show"],typing.Optional[bool]),
        save=type_formater(kwargs["save"],typing.Union[str, bool, str]),
        ax=type_formater(kwargs["ax"],typing.Optional[scanpy.plotting._utils._AxesSubplot]),
        return_fig=True,
        vmin=type_formater(kwargs["vmin"],typing.Optional[float]),
        vmax=type_formater(kwargs["vmax"],typing.Optional[float]),
        vcenter=type_formater(kwargs["vcenter"],typing.Optional[float]),
        norm=type_formater(kwargs["norm"],typing.Optional[matplotlib.colors.Normalize]),
    )

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

config.methods_plot["matrixplot"] = dict(
    
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

    function = matrixplot_plot
)
