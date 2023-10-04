import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
# from dash import dcc
from dash import dcc
from dash import html
import plotly.graph_objs as go
import plotly.figure_factory as ff
from plotly.subplots import make_subplots
import dash
import scrublet
from scipy.stats import mode
from scipy.spatial.distance import pdist, squareform
import dash_ag_grid as dag

from ..functions import *

from app import app

from .. import config

def args_differential_expression():

    options = node_names(exclude_downstream_from_node=config.selected) 
    options_layer = [i for i in config.adata.layers.keys()] 

    return [
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Use the indicated neighbors.",
            "value":None,
            "clearable":False,
            "options":options
        },
        {
            "input":"BooleanSwitch",
            "name":"use_raw",
            "description":"Use raw attribute of adata if present.",
            "value":True,
        },
        {
            "input":"Dropdown",
            "name":"layer",
            "description":"Key from adata.layers whose value will be used to perform tests on.",
            "value":None,
            "clearable":True,
            "options":options_layer,
            "summary":True
        },
        {
            "input":"Input",
            "name":"n_genes",
            "description":"The number of genes that appear in the returned tables. Defaults to all genes.",
            "value":None,
            "type":"number"
        },
        {
            "input":"Dropdown",
            "name":"method",
            "description":"The default method is 't-test', 't-test_overestim_var' overestimates variance of each group, 'wilcoxon' uses Wilcoxon rank-sum, 'logreg' uses logistic regression. See [Ntranos18], here and here, for why this is meaningful.",
            "value":"t-test",
            "clearable":False,
            "options":["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
            "summary":True
        },
        {
            "input":"Dropdown",
            "name":"corr_method",
            "description":"p-value correction method. Used only for 't-test', 't-test_overestim_var', and 'wilcoxon'.",
            "value":"benjamini-hochberg",
            "clearable":False,
            "options":["benjamini-hochberg", "bonferroni"],
            "summary":True
        },
        {
            "input":"BooleanSwitch",
            "name":"tie_correct",
            "description":"Use tie correction for 'wilcoxon' scores. Used only for 'wilcoxon'.",
            "value":False,
        },
        {
            "input":"BooleanSwitch",
            "name":"pts",
            "description":"Compute the fraction of cells expressing the genes.",
            "value":False,
        },
    ]

def f_differential_expression(name_analysis, kwargs):

    if "log1p" not in config.adata.uns.keys(): #Solve in case there is a problem with it
        config.adata.uns["log1p"] = {"base":10}
    elif "base" not in config.adata.uns["log1p"].keys():
        config.adata.uns["log1p"] = {"base":10}

    sc.tl.rank_genes_groups(
                config.adata,
                groupby=kwargs["input"],
                use_raw=kwargs["use_raw"],
                layer=kwargs["layer"],
                n_genes=kwargs["n_genes"],
                method=kwargs["method"],
                corr_method=kwargs["corr_method"],
                tie_correct=kwargs["tie_correct"],
                pts=kwargs["pts"],
                key_added=name_analysis
                )
    
    pos = get_node_pos(name_analysis)
    config.graph[pos]["data"]["plotting"] = {"n_genes":2,"style":"table","cluster":config.adata.uns[config.selected]["scores"].dtype.names[0]}

def rm_differential_expression(name_analysis):

    del config.adata.uns[name_analysis]

    return

def rename_differential_expression(name_analysis, name_new_analysis):

    config.adata.uns[name_new_analysis] = config.adata.uns[name_analysis]
    del config.adata.uns[name_analysis]

    return

def plot_differential_expression(name_analysis):
    
    node = get_node(name_analysis)

    if not node["data"]["computed"]:
        return []

    l = []

    if node["data"]["plotting"]["style"] == "heatmap":
        l = heatmap(node["data"]["plotting"]["n_genes"])
    elif node["data"]["plotting"]["style"] == "clustermap":
        l = clustermap(node["data"]["plotting"]["n_genes"])
    elif node["data"]["plotting"]["style"] == "table":
        l = table()

    return l

@app.callback(
    dash.Output("analysis_use_raw","on", allow_duplicate=True),
    dash.Input("analysis_input","value"),
    prevent_initial_call=True
)
def change_raw(v):

    if config.adata.raw:
        return True
    else:
        return False
    
def clustermap(n_genes):

    name_analysis = config.selected

    clusters = config.adata.uns[name_analysis]["scores"].dtype.names
    labels_x = []
    for cluster in clusters:
        order = np.argsort(-np.abs(config.adata.uns[name_analysis]["scores"][cluster]))
        labels_x = np.append(labels_x, config.adata.uns[name_analysis]["names"][cluster][order[:n_genes]])

    labels_y = clusters

    data_array = np.zeros([n_genes*len(clusters),len(clusters)])
    color = []
    for i,cluster in enumerate(clusters):
        for j,gene in enumerate(labels_x):
            g = gene == config.adata.uns[name_analysis]["names"][cluster]
            data_array[j,i] = config.adata.uns[name_analysis]["scores"][cluster][g][0]
            color.append(config.adata.uns[name_analysis]["scores"][cluster][g][0])

    ratio = int(data_array.shape[0]/data_array.shape[1])

    fig = make_subplots(rows=2, cols=2, 
                        column_widths=[1-0.2/ratio, 0.2/ratio], row_heights=[0.2, 0.8],
                        horizontal_spacing=0, vertical_spacing=0,
                        shared_xaxes=True, shared_yaxes=True,
                        )

    # Initialize figure by creating upper dendrogram
    p1 = list(range(0, data_array.shape[0], ratio))
    data = data_array[p1,:]
    for i in range(1,int(data_array.shape[0]/data_array.shape[1])):
        p1 = list(range(0, data_array.shape[0], ratio))
        data = np.append(data,data_array[p1,:],axis=1)

    dendro_top = ff.create_dendrogram(data, orientation='bottom')#, labels=labels_y)
    mm = np.inf
    MM = -np.inf
    for i in range(len(dendro_top['data'])):
        m = min(dendro_top['data'][i]['x'])
        M = max(dendro_top['data'][i]['x'])
        mm = min(m,mm)
        MM = max(M,MM)
    for i in range(len(dendro_top['data'])):
        dendro_top['data'][i]['yaxis'] = 'y2'
        dendro_top['data'][i]['x'] -= mm
        dendro_top['data'][i]['x'] *= (len(labels_x)-ratio)/(MM-mm)
        dendro_top['data'][i]['x'] += ratio/4
    
    dendro_leaves_x =[]
    for i in np.array(dendro_top['layout']['xaxis']['ticktext'],int)*ratio:
        for j in range(ratio):
            dendro_leaves_x.append(i+j)

    # Add Side Dendrogram Data to Figure
    for data in dendro_top['data']:
        fig.add_trace(data, row=1, col=1)
    
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(showticklabels=False, row=1, col=1)

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array.transpose(), orientation='left')
    mm = np.inf
    MM = -np.inf
    for i in range(len(dendro_side['data'])):
        m = min(dendro_side['data'][i]['y'])
        M = max(dendro_side['data'][i]['y'])
        mm = min(m,mm)
        MM = max(M,MM)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
        dendro_side['data'][i]['y'] -= mm
        dendro_side['data'][i]['y'] *= (len(labels_y)-1)/(MM-mm)

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data, row=2, col=2)

    fig.update_xaxes(showticklabels=False, row=2, col=2)
    fig.update_yaxes(showticklabels=False, row=2, col=2)

    # Create Heatmap
    dendro_leaves_y = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves_y = list(map(int, dendro_leaves_y))
    heat_data = data_array.copy()
    # heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves_x,:]
    heat_data = heat_data[:,dendro_leaves_y]
    heat_data = heat_data.transpose()

    X,Y = np.meshgrid(np.arange(0,len(labels_y),1),np.arange(0,len(labels_x),1))
    X = X[dendro_leaves_x,:]
    X = X[:,dendro_leaves_y]
    Y = Y[dendro_leaves_x,:]
    Y = Y[:,dendro_leaves_y]

    heatmap = [
        go.Scatter(
            x = Y.reshape(-1),
            y = X.reshape(-1),
            mode = "markers",
            marker = {
                "color":heat_data.reshape(-1),
                "size":np.abs(heat_data.reshape(-1)),
            },
            text=list(map(str,heat_data.reshape(-1)))
        )
    ]

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data, row=2, col=1)

    fig.update_xaxes(tickvals=np.arange(0,len(labels_x),1),
                     ticktext=np.array(labels_x)[dendro_leaves_x],
                     row=2, col=1)
    fig.update_yaxes(tickvals=np.arange(0,len(labels_y),1),
                     ticktext=np.array(labels_y)[dendro_leaves_y],
                     row=2, col=1)

    # Edit Layout
    fig.update_layout({'width':700*ratio, 'height':800,
                            'showlegend':False, 'hovermode': 'closest',
                            })

    return [
            dcc.Dropdown(
                id = "differential_expression_plot_style",
                value=get_node(config.selected)["data"]["plotting"]["style"],
                options=["heatmap","clustermap","table"],
                clearable=False
            ),
            dcc.Dropdown(
                id = "differential_expression_plot_n_genes",
                value=get_node(config.selected)["data"]["plotting"]["n_genes"],
                options=[i for i in range(1,10)],
                clearable=False
            ),
            dbc.Col(),
            dbc.Col(
                dcc.Graph(figure=fig)
            ),
            dbc.Col(),
    ]

def heatmap(n_genes):

    name_analysis = config.selected

    clusters = config.adata.uns[name_analysis]["scores"].dtype.names
    labels_x = []
    for cluster in clusters:
        order = np.argsort(-np.abs(config.adata.uns[name_analysis]["scores"][cluster]))
        labels_x = np.append(labels_x, config.adata.uns[name_analysis]["names"][cluster][order[:n_genes]])

    labels_y = clusters

    data_array = np.zeros([n_genes*len(clusters),len(clusters)])
    color = []
    for i,cluster in enumerate(clusters):
        for j,gene in enumerate(labels_x):
            g = gene == config.adata.uns[name_analysis]["names"][cluster]
            data_array[j,i] = config.adata.uns[name_analysis]["scores"][cluster][g][0]
            color.append(config.adata.uns[name_analysis]["scores"][cluster][g][0])

    ratio = int(data_array.shape[0]/data_array.shape[1])

    fig = make_subplots(rows=2, cols=2, 
                        column_widths=[1-0.2/ratio, 0.2/ratio], row_heights=[0.2, 0.8],
                        horizontal_spacing=0, vertical_spacing=0,
                        shared_xaxes=True, shared_yaxes=True,
                        )

    # Initialize figure by creating upper dendrogram
    p1 = list(range(0, data_array.shape[0], ratio))
    data = data_array[p1,:]
    for i in range(1,int(data_array.shape[0]/data_array.shape[1])):
        p1 = list(range(0, data_array.shape[0], ratio))
        data = np.append(data,data_array[p1,:],axis=1)

    dendro_top = ff.create_dendrogram(data, orientation='bottom')#, labels=labels_y)
    mm = np.inf
    MM = -np.inf
    for i in range(len(dendro_top['data'])):
        m = min(dendro_top['data'][i]['x'])
        M = max(dendro_top['data'][i]['x'])
        mm = min(m,mm)
        MM = max(M,MM)
    for i in range(len(dendro_top['data'])):
        dendro_top['data'][i]['yaxis'] = 'y2'
        dendro_top['data'][i]['x'] -= mm
        dendro_top['data'][i]['x'] *= (len(labels_x)-ratio)/(MM-mm)
        dendro_top['data'][i]['x'] += ratio/4
    
    dendro_leaves_x =[]
    for i in np.array(dendro_top['layout']['xaxis']['ticktext'],int)*ratio:
        for j in range(ratio):
            dendro_leaves_x.append(i+j)

    # Add Side Dendrogram Data to Figure
    for data in dendro_top['data']:
        fig.add_trace(data, row=1, col=1)
    
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(showticklabels=False, row=1, col=1)

    # Create Side Dendrogram
    dendro_side = ff.create_dendrogram(data_array.transpose(), orientation='left')
    mm = np.inf
    MM = -np.inf
    for i in range(len(dendro_side['data'])):
        m = min(dendro_side['data'][i]['y'])
        M = max(dendro_side['data'][i]['y'])
        mm = min(m,mm)
        MM = max(M,MM)
    for i in range(len(dendro_side['data'])):
        dendro_side['data'][i]['xaxis'] = 'x2'
        dendro_side['data'][i]['y'] -= mm
        dendro_side['data'][i]['y'] *= (len(labels_y)-1)/(MM-mm)

    # Add Side Dendrogram Data to Figure
    for data in dendro_side['data']:
        fig.add_trace(data, row=2, col=2)

    fig.update_xaxes(showticklabels=False, row=2, col=2)
    fig.update_yaxes(showticklabels=False, row=2, col=2)

    # Create Heatmap
    dendro_leaves_y = dendro_side['layout']['yaxis']['ticktext']
    dendro_leaves_y = list(map(int, dendro_leaves_y))
    heat_data = data_array.copy()
    # heat_data = squareform(data_dist)
    heat_data = heat_data[dendro_leaves_x,:]
    heat_data = heat_data[:,dendro_leaves_y]
    heat_data = heat_data.transpose()

    heatmap = [
        go.Heatmap(
            x = dendro_leaves_x,
            y = dendro_leaves_y,
            z = heat_data,
            colorscale = 'Blues'
        )
    ]

    heatmap[0]['x'] = np.arange(0,len(labels_x),1)
    heatmap[0]['y'] = np.arange(0,len(labels_y),1)

    # Add Heatmap Data to Figure
    for data in heatmap:
        fig.add_trace(data, row=2, col=1)

    fig.update_xaxes(tickvals=np.arange(0,len(labels_x),1),
                     ticktext=np.array(labels_x)[dendro_leaves_x],
                     row=2, col=1)
    fig.update_yaxes(tickvals=np.arange(0,len(labels_y),1),
                     ticktext=np.array(labels_y)[dendro_leaves_y],
                     row=2, col=1)

    # Edit Layout
    fig.update_layout({'width':700*ratio, 'height':800,
                            'showlegend':False, 'hovermode': 'closest',
                            })

    return [
            dcc.Dropdown(
                id = "differential_expression_plot_style",
                value=get_node(config.selected)["data"]["plotting"]["style"],
                options=["heatmap","clustermap","table"],
                clearable=False
            ),
            dcc.Dropdown(
                id = "differential_expression_plot_n_genes",
                value=get_node(config.selected)["data"]["plotting"]["n_genes"],
                options=[i for i in range(1,10)],
                clearable=False
            ),
            dbc.Col(),
            dbc.Col(
                dcc.Graph(figure=fig)
            ),
            dbc.Col(),
    ]

def table():

    l = ["names","scores","pvals","pvals_adj"]

    df = pd.DataFrame(columns = l,
                     data={i:config.adata.uns[config.selected][i][get_node(config.selected)["data"]["plotting"]["cluster"]] for i in l}
    )

    fig = dag.AgGrid(
        id="row-sorting-simple",
        rowData=df.to_dict("records"),
        columnDefs=[{"field":i} for i in l],
        defaultColDef={"resizable": True, "sortable": True, "filter": False},
        columnSize="sizeToFit",
    )

    return [
            dcc.Dropdown(
                id = "differential_expression_plot_style",
                value=get_node(config.selected)["data"]["plotting"]["style"],
                options=["heatmap","clustermap","table"],
                clearable=False
            ),
            dcc.Dropdown(
                id = "differential_expression_cluster",
                value=get_node(config.selected)["data"]["plotting"]["cluster"],
                options=config.adata.uns[config.selected]["scores"].dtype.names,
                clearable=False
            ),
            dbc.Col(
                fig
            ),
    ]

@app.callback(
    dash.Output("analysis_plot","children", allow_duplicate=True),
    dash.Input("differential_expression_plot_style", "value"),
    prevent_initial_call = True
)
def change_style(val):

    prevent_race("differential_expression")

    pos = get_node_pos(config.selected)
    config.graph[pos]["data"]["plotting"]["style"] = val

    return plot_differential_expression(config.selected)

@app.callback(
    dash.Output("analysis_plot","children", allow_duplicate=True),
    dash.Input("differential_expression_plot_n_genes", "value"),
    prevent_initial_call = True
)
def change_style(val):

    prevent_race("differential_expression")

    pos = get_node_pos(config.selected)
    config.graph[pos]["data"]["plotting"]["n_genes"] = val

    return plot_differential_expression(config.selected)

@app.callback(
    dash.Output("analysis_plot","children", allow_duplicate=True),
    dash.Input("differential_expression_cluster", "value"),
    prevent_initial_call = True
)
def change_style(val):

    prevent_race("differential_expression")

    pos = get_node_pos(config.selected)
    config.graph[pos]["data"]["plotting"]["cluster"] = val

    return plot_differential_expression(config.selected)