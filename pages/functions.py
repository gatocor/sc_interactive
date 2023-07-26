import scanpy as sc
import numpy as np
import pandas as pd
import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_daq as daq

def f_qc(adata):
    if "X_raw" not in adata.obsm.keys():
        adata.obsm["X_raw"] = adata.X.copy()
    if "qc_total_counts" not in adata.obs.columns.values:
        # sc.pp.calculate_qc_metrics(adata,percent_top=(3,),inplace=True)#np.round(np.int,np.linspace(1,adata.shape[1],5)))
        adata.obs["qc_total_counts"] = np.array(adata.obsm["X_raw"].sum(axis=1)).reshape(-1)
        adata.obs["qc_n_genes_by_counts"] = np.array((adata.obsm["X_raw"] > 0).sum(axis=1)).reshape(-1)
    if "qc" not in adata.uns.keys():
        adata.uns["qc"] = {
            "qc_total_counts":{"Minimum threshold":adata.obs["qc_total_counts"].min(), "Maximum threshold":adata.obs["qc_total_counts"].max()},
            "qc_n_genes_by_counts":{"Minimum threshold":adata.obs["qc_n_genes_by_counts"].min(), "Maximum threshold":adata.obs["qc_n_genes_by_counts"].max()},
        }
    if "gene_lists" not in adata.uns.keys():
        adata.uns["gene_lists"] = {}

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
                    # adata.obs["qc_"+i["Pattern"]] = np.array(adata.obsm["X_raw"][:,pp].sum(axis=1)).reshape(-1)/adata.obs["qc_total_counts"]
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

def f_qc_table_threshold(adata):
    
    l = [{"Control measure":i,
          "Minimum Threshold":str(adata.uns["qc"][i]["Minimum threshold"]),
          "Maximum Threshold":str(adata.uns["qc"][i]["Maximum threshold"]),
          } for i in ["qc_total_counts","qc_n_genes_by_counts"]]

    l_pattern = [{"Control measure":i,
          "Minimum Threshold":str(adata.uns["qc"]["patterns"][i]["Minimum threshold"]),
          "Maximum Threshold":str(adata.uns["qc"]["patterns"][i]["Maximum threshold"]),
          } for i in adata.uns["qc"]["patterns"]]

    return l+l_pattern

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
