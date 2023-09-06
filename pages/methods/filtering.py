import numpy as np
import scanpy as sc
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import dash
import dash_table

from app import app

from .. import config

def filtering_args(adata):

    options = []
    optionsthreshold = []
    if "__interactive__" in adata.uns.keys():
        options = [i for i,j in adata.uns["__interactive__"].items() if j["type"] == "QC"] 
        optionsthreshold = ["total_counts","n_genes_by_counts"]
    
    return [
        "Cell filtering",
        {
            "input":"Dropdown",
            "name":"input",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options":options
        },
        {
            "input":"ThresholdTable",
            "name":"thresholds",
            "description":"Metrics of quality control to compute",
            "value":[],
            "options":optionsthreshold
        },
    ]

def f_filtering(adata, name_analysis, **kwargs):
        
    return

def make_filtering_plots(adata, name_analysis):

    l = []
    dic = adata.uns["__interactive__"][name_analysis]["params"]["thresholds"] 

    # l += [dbc.Row(
    #     dash_table.DataTable(
    #             id="bins",
    #             columns=[
    #                 {"name": i, "id": i, "deletable": False, "editable": True if ("#bins" in i) else False} for i in ["Measure","#bins"]
    #             ],
    #             data=[{"Measure":i["Measure"],"#bins":30} for i in dic],
    #             editable=False, 
    #             row_deletable=False,
    #             style_table={'overflowY': 'auto', 'overflowX': 'auto'},
    #             style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
    #     ),
    # )]
    
    for metric in dic:

        print(name_analysis)

        var_selected_data = metric["Measure"]
            
        hist_values, hist_bins = np.histogram(adata.obs[var_selected_data].values, bins=int(np.round(int(metric["#bins"]))))
        tallest_bin_height = np.max(hist_values)
        var1_vertical_line_min = go.Scatter(
            x=[metric["Min"], metric["Min"]],
            y=[0, tallest_bin_height],
            mode='lines',
            name='Min threshold',
            line=dict(color='red', dash='dash')
        )
        var1_vertical_line_max = go.Scatter(
            x=[metric["Max"], metric["Max"]],
            y=[0, tallest_bin_height],
            mode='lines',
            name='Max threshold',
            line=dict(color='green', dash='dash')
        )

        l += [
                dbc.Row(
                    [dcc.Graph(id="Histogram",
                            figure={
                                    "data":[
                                        go.Histogram(
                                            x=adata.obs[var_selected_data].values,
                                            nbinsx=int(np.round(int(metric["#bins"]))),
                                            name='Histogram',
                                            marker=dict(color='blue'),
                                            opacity=0.7
                                        ),
                                        var1_vertical_line_min,
                                        var1_vertical_line_max
                                    ],
                                    "layout":{
                                            'title': f'Histogram of {var_selected_data}',
                                            'xaxis': {'title': var_selected_data},
                                            'yaxis': {'title': 'Count'},
                                            'barmode': 'overlay',
                                            'width':1500,
                                            'height':400,
                                    }
                                }
                    )],
                    justify="center",
                    style={'width': '90%', 'margin': 'auto'}
                )
            ]
        
    return l

# @app.callback(
#     dash.Output("filtering_measure","data"),
#     dash.Input("filtering_measure_button","n_clicks"),
#     dash.State("filtering_measure_dropdown","value"),
#     dash.State("filtering_measure","data"),
#     prevent_initial_callback=True
# )
# def add_metric(_,value,data):

#     values = [i["Measure"] for i in data]
#     if value != None and value not in values:
#         data.append({"Measure":value})

#     return data

@app.callback(
    dash.Output("filtering_thresholds","data"),
    dash.Input("filtering_input","value"),
    prevent_initial_callback=True
)
def show_table(value):

    d = []
    if value != None:
        l = [i for i,j in config.adata.uns["__interactive__"].items() if j["type"]=="QC"]
        d = []
        for i in l:
            for j in config.adata.uns["__interactive__"][i]["params"]["measure"][1:-1].split(","):
                d.append({"Measure":j,"Min":0,"Max":0,"#bins":30})

    if "params" not in config.adata.uns["__interactive__"][config.name_analysis].keys():
        config.adata.uns["__interactive__"][config.name_analysis]["params"] = {}

    config.adata.uns["__interactive__"][config.name_analysis]["params"]["thresholds"] = d

    return d

@app.callback(
    dash.Output("filtering_plots2","children",allow_duplicate=True),
    dash.Input("filtering_thresholds","data"),
    prevent_initial_call=True
)
def update_table(value):

    config.adata.uns["__interactive__"][config.name_analysis]["params"]["thresholds"] = value

    return make_filtering_plots2(config.adata, config.name_analysis)