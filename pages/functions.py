import scanpy as sc
import numpy as np
import pandas as pd
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq
import dash_table

def f_qc_base(adata):
    if "X_raw" not in adata.obsm.keys():
        adata.obsm["X_raw"] = adata.X.copy()
    if "total_counts" not in adata.obs.columns.values:
        adata.obs["total_counts"] = np.array(adata.obsm["X_raw"].sum(axis=1)).reshape(-1)
    if "n_genes_by_counts" not in adata.obs.columns.values:
        adata.obs["n_genes_by_counts"] = np.array((adata.obsm["X_raw"] > 0).sum(axis=1)).reshape(-1)
    if "qc" not in adata.uns.keys():
        adata.uns["qc"] = {
            "total_counts":{"Minimum threshold":adata.obs["total_counts"].min(), "Maximum threshold":adata.obs["total_counts"].max()},
            "n_genes_by_counts":{"Minimum threshold":adata.obs["n_genes_by_counts"].min(), "Maximum threshold":adata.obs["n_genes_by_counts"].max()},
        }
    if "total_counts" not in adata.uns["qc"].keys():
        adata.uns["qc"]["total_counts"] = {"Minimum threshold":adata.obs["total_counts"].min(), "Maximum threshold":adata.obs["total_counts"].max()}
    if "n_genes_by_counts" not in adata.uns["qc"].keys():
        adata.uns["qc"]["n_genes_by_counts"] = {"Minimum threshold":adata.obs["n_genes_by_counts"].min(), "Maximum threshold":adata.obs["n_genes_by_counts"].max()}
    if "gene_lists" not in adata.uns.keys():
        adata.uns["gene_lists"] = {}

def f_qc(adata,metrics,data):

    m = [i["Name"] for i in metrics]
    l = list(adata.uns["qc"].keys())

    if "total_counts" not in adata.obs.columns.values:
        adata.obs["total_counts"] = np.array(adata.obsm["X_raw"].sum(axis=1)).reshape(-1)
    if "n_genes_by_counts" not in adata.obs.columns.values:
        adata.obs["n_genes_by_counts"] = np.array((adata.obsm["X_raw"] > 0).sum(axis=1)).reshape(-1)
    if "total_counts" not in m:
        adata.uns["qc"]["total_counts"] = {"Minimum threshold":adata.obs["total_counts"].min(), "Maximum threshold":adata.obs["total_counts"].max()}
    if "n_genes_by_counts" not in m:
        adata.uns["qc"]["n_genes_by_counts"] = {"Minimum threshold":adata.obs["n_genes_by_counts"].min(), "Maximum threshold":adata.obs["n_genes_by_counts"].max()}
    
    for i in l:
        if (i not in ["total_counts","n_genes_by_counts"]) and (i not in m):
            del adata.uns["qc"][i]

    for i in m:
        if (i not in ["total_counts","n_genes_by_counts"]) and (i not in l):
            adata.obs[i] = np.array(adata.obsm["X_raw"][:,adata.var[i].values].sum(axis=1)).reshape(-1)/adata.obs["total_counts"]
            adata.uns["qc"][i] = {"Minimum threshold":adata.obs[i].min(), "Maximum threshold":adata.obs[i].max()}
        else:
            for j in data:
                if j["Control measure"] == i:
                    adata.uns["qc"][i]["Minimum threshold"] = j["Minimum threshold"]
                    adata.uns["qc"][i]["Maximum threshold"] = j["Maximum threshold"]

    return

def f_options(adata,motive):
    return [i for i in adata.obs.columns.values if i.startswith(motive)]

def f_update_patterns(adata, pattern):

    p = adata.uns["gene_lists"]
    for i in pattern:
        if i["Concept"] == '' and i["Pattern"] == '' and i["Genes"] == '':
            adata.uns["gene_lists"][""] = {"Pattern":"", "Genes":[]}
        else:
            if i["Concept"] != '':
                if i["Pattern"] != '' and i["Pattern"] != None:
                    pp = [j.startswith(i["Pattern"]) for j in adata.var["Gene"].values]
                    # adata.obs[""+i["Pattern"]] = np.array(adata.obsm["X_raw"][:,pp].sum(axis=1)).reshape(-1)/adata.obs["total_counts"]
                    adata.uns["gene_lists"][i["Concept"]] = {"Pattern":i["Pattern"], "Genes":list(adata.var["Gene"].values[pp])}
                elif "[" not in i["Genes"]:
                    adata.uns["gene_lists"][i["Concept"]] = {"Pattern":i["Pattern"], "Genes":list([j for j in i["Genes"].split(" ")])}
                else:
                    adata.uns["gene_lists"][i["Concept"]] = {"Pattern":i["Pattern"], "Genes":adata.uns["gene_lists"][i["Concept"]["Genes"]]}
    
    for i in [j for j in p.keys()]:
        if i not in [j["Concept"] for j in pattern]:
            del p[i]
            adata.var.drop(i,axis=1,inplace=True)
        elif i != '':
            adata.var[i] = [j in adata.uns["gene_lists"][i]["Genes"] for j in adata.var["Gene"].values]

    return

def f_qc_table_pattern(adata):

    return [{"Concept":i, "Pattern":str(adata.uns["gene_lists"][i]["Pattern"]),"Genes":str(adata.uns["gene_lists"][i]["Genes"])} for i in adata.uns["gene_lists"]]

def f_qc_table_metrics(adata):
    
    l = [{"Name":i} for i in adata.uns["qc"]]

    return l

def f_qc_table_threshold(adata):
    
    l = [{"Control measure":i,
          "Minimum threshold":str(adata.uns["qc"][i]["Minimum threshold"]),
          "Maximum threshold":str(adata.uns["qc"][i]["Maximum threshold"]),
          } for i in adata.uns["qc"]]

    return l

def f_qc_summary(adata):
    statistics = {}
    # statistics["Type"] = np.array(adata.uns["qc"].keys())
    # statistics["Threshold min"] = [adata.uns["qc"]["patterns"][i]["Minimum threshold"] for i in adata.uns["qc"]["patterns"]]
    # statistics["Threshold max"] = [adata.uns["qc"]["patterns"][i]["Maximum threshold"] for i in adata.uns["qc"]["patterns"]]
    # statistics["Cells rejected"] = np.array([np.sum((adata.obs[i].values<adata.uns["qc"]["patterns"][i]["Minimum threshold"]) + (adata.obs[i].values>adata.uns["qc"]["patterns"][i]["Maximum threshold"]) > 0) for i in adata.uns["qc"].keys()])
    # statistics["Cells rejected %"] = statistics["Cells rejected"]/adata.shape[0]

    return pd.DataFrame(statistics)

def make_arguments(id, arglist, title="Arguments", add_execution_button=True):

    l = []

    l.append(html.H1(title))
    for i,arg in enumerate(arglist):
        l.append(
            dbc.Tooltip(
                    arg["description"],
                    target=id+str(i),
                    placement="bottom",
            )
        )
        lab = html.Label(arg["name"],id=id+str(i))
        if arg["input"] == "Input":
            input = dbc.Input(id=id+"_"+str(i),value=arg["value"],type=arg["type"])
        elif arg["input"] == "Dropdown":
            input = dcc.Dropdown(
                        id=id+"_"+str(i),
                        options=arg["options"],
                        value=arg["value"],
                        # placeholder="Select a column",
                        clearable=False
                    )
        elif arg["input"] == "BooleanSwitch":
            input = daq.BooleanSwitch(id=id+"_"+str(i), on=arg["value"])
        else:
            print(arg["input"]," input does not exist.")

        l.append(
            dbc.Row(
                [
                    dbc.Col(
                        lab,
                        # width=text_width
                    ),
                    dbc.Col(
                        input,
                        # width=input_width
                    )
                ],
            )
        )

    if add_execution_button:
        l.append(
            dbc.Row(
                dbc.Button("Execute",id="button_"+id,                           )
            )
        )

    return l

def json_serializable(uns):

    d = {}

    for i in uns:
        if type(uns[i]) == dict:
            d[i] = json_serializable(uns[i])
        elif type(uns[i]) in [int, float, str, list, np.int_, np.float_]:
            d[i] = uns[i]
        else:
            d[i] = str(type(uns[i]))

    return d

def make_qc_per_condition(adata):

    print(adata.uns["qc"]["total_counts"].keys())
    if "per_condition" in adata.uns["qc"]["total_counts"].keys():

        for var_selected_data in adata.uns["qc"].keys():          
            del adata.uns["qc"][var_selected_data]["per_condition"]

        return [dbc.Button(id='qc_per_condition-button', n_clicks=0, children="Add per conditon analysis",
                        size="lg",
                        style={
                            "background-color":"gray",
                            # 'width': '680px', 
                        }      
                        )]
    else:
        
        for var_selected_data in adata.uns["qc"].keys():          
            adata.uns["qc"][var_selected_data]["per_condition"] = {}

        return  [
                dbc.Row(
                    [
                        dbc.Col(
                            html.H1("Per condition Analysis")
                        ),
                        dbc.Col(
                            dbc.Button(id='qc_per_condition-button', n_clicks=1, children="Remove per condition analysis",
                                    size="lg",
                                    style={
                                        "background-color":"#343A40",
                                        'width': '280px', 
                                    }      
                                    ),
                            width=2
                        ),
                    ],
                    style={"margin-bottom":"1cm"}
                ),
                dbc.Row(
                        dash_table.DataTable(
                            id='table_qc_per_condition_metrics',
                            columns=[
                                {"name": i, "id": i, "deletable": False, "editable": False} for i in ["Name"]
                            ],
                            data=[],
                            editable=False,
                            row_deletable=True,
                            style_table={'overflowY': 'auto', 'overflowX': 'auto'},
                            style_cell={'textAlign': 'left', 'whiteSpace': 'normal', 'height': 'auto'}
                        ),
                        # style={"margin-bottom":"1cm"}        
                    ),
                dbc.Row(
                    [
                        dbc.Col(dcc.Dropdown(
                            id = "dropdown_add_per_condition_metrics",
                            value = None,
                            options = [str(i) for i in adata.obs.columns.values if (adata.obs.dtypes[i] in ["category" ,object, str, int])]
                        )),
                        dbc.Col(dbc.Button("Add condition", id="add_qc_per_condition-button")),
                    ],
                    justify="Left",
                    className="mb-4",
                ),
                dbc.Row(
                    id="per_condition_plot",
                    children=dcc.Graph(),
                )
            ]