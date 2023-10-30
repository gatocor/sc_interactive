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

from app import app

def plot_center(fig):
    return [
        dbc.Col(),
        dbc.Col(fig),
        dbc.Col()
    ]

def plot_violin(x=None, y=None, bandwidth=1):

    return go.Violin(
                        x=x,
                        y=y,
            )

def plot_scatter(x=None, y=None, color=None, marker_size=1):

    if color == None:
        marker = {}
    else:
        marker = {"color":qualitative_colors(color),"size":marker_size}

    return go.Scatter(
                        x=x,
                        y=y,
                        mode="markers",
                        marker=marker
            )

def make_Graph(data=[],title="",x_label="",y_label="",style={"width": "50vw", "height": "50vh"}, center=True):
     
    graph =  dcc.Graph(
        figure={
                    "data":data,
                    "layout":{
                            "title": title,
                            "xaxis": {"title": x_label},
                            "yaxis": {"title": y_label},
                    }
        },
        style=style
    )

    if center:

        graph = [
            dbc.Col(),
            dbc.Col(graph),
            dbc.Col()
        ]

    return graph

def plot_dendogram(data_array, groups=None, orientation='bottom'):

    if orientation in ["bottom", "top"]:
        v = "x"
        v2 = "y"
    else:
        v = "y"
        v2 = "x"
        data_array = data_array.transpose()

    ratio = 1
    if groups == None:
        ratio = int(data_array.shape[0]/data_array.shape[1])

    p1 = list(range(0, data_array.shape[0], ratio))
    data = data_array[p1,:]
    for i in range(1,ratio):
        p1 = list(range(0, data_array.shape[0], ratio))
        data = np.append(data,data_array[p1,:],axis=1)

    dendro_top = ff.create_dendrogram(data, orientation=orientation)
    mm = np.inf
    MM = -np.inf
    for i in range(len(dendro_top['data'])):
        m = min(dendro_top['data'][i][v])
        M = max(dendro_top['data'][i][v])
        mm = min(m,mm)
        MM = max(M,MM)
    for i in range(len(dendro_top['data'])):
        dendro_top['data'][i][v2+'axis'] = v2+'2'
        dendro_top['data'][i][v] -= mm
        dendro_top['data'][i][v] *= (data_array.shape[0]-ratio)/(MM-mm)
        dendro_top['data'][i][v] += ratio/4
    
    dendro_leaves_x =[]
    for i in np.array(dendro_top['layout'][v+'axis']['ticktext'],int)*ratio:
        for j in range(ratio):
            dendro_leaves_x.append(i+j)

    return dendro_top, dendro_leaves_x

def plot_scattermap(data_array, xorder=None, yorder=None):

    if xorder == None:
        xorder = np.arange(0,data_array.shape[0],1)

    if yorder == None:
        yorder = np.arange(0,data_array.shape[1],1)

    data_array = data_array[xorder,:]
    data_array = data_array[:,yorder]

    X,Y = np.meshgrid(np.arange(0,len(yorder),1),np.arange(0,len(xorder),1))
    X = X[xorder,:]
    X = X[:,yorder]
    Y = Y[xorder,:]
    Y = Y[:,yorder]

    data_array = data_array.transpose()
    X = X.transpose()
    Y = Y.transpose()

    heatmap = go.Scatter(
            x = Y.reshape(-1),
            y = X.reshape(-1),
            mode = "markers",
            marker = {
                "color":data_array.reshape(-1),
                "size":np.abs(data_array.reshape(-1)),
            },
            text=list(map(str,data_array.reshape(-1)))
        )
    
    return heatmap

def plot_heatmap(data_array, xorder=None, yorder=None):

    return go.Heatmap(
            z = data_array.transpose()
        )

def plot_clustermap(data_array, labels_x, labels_y, style="heatmap"):

    fig = make_subplots(rows=2, cols=2, 
                        column_widths=[.8, .2], row_heights=[0.2, 0.8],
                        horizontal_spacing=0, vertical_spacing=0,
                        shared_xaxes=True, shared_yaxes=True,
                        )

    # Initialize figure by creating dendrograms
    dendrox, xorder = plot_dendogram(data_array, orientation='bottom')
    for data in dendrox['data']:
        fig.add_trace(data, row=1, col=1)
    fig.update_xaxes(showticklabels=False, row=1, col=1)
    fig.update_yaxes(showticklabels=False, row=1, col=1)

    dendroy, yorder = plot_dendogram(data_array, groups = 1, orientation="left")
    # dendroy, yorder = dendogram_side(data_array)
    for data in dendroy['data']:
        fig.add_trace(data, row=2, col=2)
    fig.update_xaxes(showticklabels=False, row=2, col=2)
    fig.update_yaxes(showticklabels=False, row=2, col=2)

    # Create Heatmap
    if style == "heatmap":
        scatter = plot_heatmap(data_array, xorder, yorder)
    elif style == "scattermap":
        scatter = plot_scattermap(data_array, xorder, yorder)

    # Add Heatmap Data to Figure
    fig.add_trace(scatter, row=2, col=1)
    fig.update_xaxes(tickvals=np.arange(0,len(labels_x),1),
                     ticktext=np.array(labels_x)[xorder],
                     row=2, col=1)
    fig.update_yaxes(tickvals=np.arange(0,len(labels_y),1),
                     ticktext=np.array(labels_y)[yorder],
                     row=2, col=1)

    # Edit Layout
    fig.update_layout({'width':50*data_array.shape[0], 'height':50*data_array.shape[1],
                            'showlegend':False, 'hovermode': 'closest',
                            })

    return [
            dbc.Col(),
            dbc.Col(
                dcc.Graph(figure=fig)
            ),
            dbc.Col(),
    ]

def plot_table(df, resizable=True,sortable=True,filter=False, editable=False, rowdeletable=False):

    if rowdeletable:

        l = []
        for i in df.to_dict("records"):
            i["del"] = "X"
            l.append(i)

        l.append({"del":"X","Type":"","Genes":"[]"})

        fig = dag.AgGrid(
            id="table",
            rowData=l,
            columnDefs=[{
                "field":"del", "cellRenderer": "Button", "lockPosition":'left', "cellRendererParams": {"className": "file-selector-button"},  "editable": False
            }]+[{"field":i} for i in df.columns.values],
            defaultColDef={"resizable": resizable, "sortable": sortable, "filter": filter, "editable": editable},
            deleteSelectedRows=True,
            dashGridOptions={"suppressRowClickSelection":True},
            columnSize="sizeToFit",
        )

    else:

        fig = dag.AgGrid(
            id="table",
            rowData=l,
            columnDefs=[{"field":i} for i in df.columns.values],
            defaultColDef={"resizable": resizable, "sortable": sortable, "filter": filter, "editable": editable},
            deleteSelectedRows=True,
            dashGridOptions={"suppressRowClickSelection":True},
            columnSize="sizeToFit",
        )

    return dbc.Col(
                fig
            )

