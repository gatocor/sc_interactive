import inspect
import scanpy as sc
import typing
import re

include = [
    'calculate_qc_metrics', 
    # 'combat', 
    # 'downsample_counts', 
    # 'filter_cells', 
    # 'filter_genes', 
    # 'filter_genes_dispersion', 
    # 'highly_variable_genes', 
    # 'log1p', 
    # 'neighbors', 
    # 'normalize_total', 
    # 'pca', 
    # 'recipe_seurat', 
    # 'recipe_weinreb17', 
    # 'recipe_zheng17', 
    # 'regress_out', 
    # 'scale', 
    # 'sqrt', 
    # 'subsample', 
    # 'dendrogram', 
    # 'diffmap', 
    # 'dpt', 
    # 'draw_graph', 
    # 'embedding_density', 
    # 'filter_rank_genes_groups', 
    # 'ingest', 
    # 'leiden', 
    # 'louvain', 
    # 'marker_gene_overlap', 
    # 'paga', 
    # 'paga_compare_paths', 
    # 'paga_degrees', 
    # 'paga_expression_entropies', 
    # 'pca', 
    # 'rank_genes_groups', 
    # 'score_genes', 
    # 'score_genes_cell_cycle', 
    # 'tsne', 
    # 'umap'         
]
models = [("pp",i) for i in dir(sc.pp) if i in include]
models += [("tl",i) for i in dir(sc.tl) if i in include]
# models += [("epp",i) for i in dir(sc.external.pp) if not i.startswith("_") and i not in ["base","utils"]]
# models += [("etl",i) for i in dir(sc.external.tl) if not i.startswith("_") and i not in ["base","utils"]]

plot = [("pl",i) for i in dir(sc.pl) if not i.startswith("_")]

print([i[1] for i in models])

# for model in models:
#     print(model)
COUNT = 0
METHOD = ""

def adapt(args, args_info):

    global COUNT, METHOD

    kargs = []
    txt_missing = ""
    txt_ignored = ""
    n_args = 0
    if args.defaults != None:
        n_args = len(args.args)-len(args.defaults)
        if n_args > 1:
            txt_missing += "\n\t\t"+str(args.args[:n_args])
        for i,d in enumerate(args.defaults):
            try:
                kargs.append((args.args[i+n_args],d,args.annotations[args.args[i+n_args]]))
            except:
                kargs.append((args.args[i+n_args],d,None))

    if args.kwonlyargs != None:
        for i,d in enumerate(args.kwonlyargs):
            kargs.append((d,args.kwonlydefaults[d],args.annotations[d]))

    executioncode = ""
    codeargs = ""
    for j,arg in enumerate(kargs):
        if arg == "self":
            None
        elif arg == "adata":
            None
        else:
            try:
                res = re.findall(arg+r'\n.*\n', args_info)[0].split("\n")[1].split("  ")[-1]
            except:
                res = ""
            res = res.replace('`','\'')
            res = res.replace('\'"','\'')
            res = res.replace('"\'','\'')
            res = arg[2]

            value = arg[1]

            codearg = f"""
        {arg[0]}=kwargs['{arg[0]}'],"""

            if arg[0] == "batch_key":
                method = f"""
    dict(
        input = "Dropdown",
        name = "{arg[0]}",
        description = "{res}",
        properties = dict(
            value = None,
            clearable = True,
            options = dict(function = "[i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in ['str','int','categorical']]")
        )
    )"""
            elif arg[0] == "layer" or arg[0] == "layers":
                method = f"""
    dict(
        input = "Dropdown",
        name = "{arg[0]}",
        description = "{res}",
        properties = dict(
            value = None,
            clearable = True,
            options = dict(function = "['X']+[i for i in config.adata.layers.keys()]")
        )
    )"""
            elif arg[0] == "init_pos" and METHOD == "umap":
                method = f"""
    dict(
        input = "Dropdown",
        name = "init_pos",
        description = "{res}",
        properties = dict(
            value = 'spectral',
            clearable = True,
            options = ['paga', 'spectral', 'random']
        )
    )"""
            elif arg[0] == "inplace" and METHOD == "calculate_qc_metrics":
                method = f"""
    dict(
        input = "BooleanSwitch",
        name = "inplace",
        description = "{res}",
        properties = dict(
            on = True,
            disabled=True
        )
    )"""
            elif arg[0] == "svd_solver" and METHOD == "recipe_weinreb17":
                method = f"""
    dict(
        input = "Dropdown",
        name = "init_pos",
        description = "{res}",
        properties = dict(
            value = 'randomized',
            clearable = True,
            options = ['randomized']
        )
    )"""
            elif arg[0] == "qc_vars" and METHOD == "calculate_qc_metrics":

                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value="",type="text")
    )"""
                codearg = f"""
        {arg[0]}=type_formater(kwargs["{arg[0]}"],typing.Union[tuple,list],typing.Union[str]),"""
                
            elif arg[0] == "percent_top" and METHOD == "calculate_qc_metrics":

                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value="(50, 100, 200, 500)",type="text"),
    )"""
                codearg = f"""
        {arg[0]}=type_formater(kwargs["{arg[0]}"],typing.Union[tuple,list],typing.Union[int,float]),"""
                
            elif arg[0] == "covariates" and METHOD == "combat":
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value="",type="text")
    )"""
            elif arg[0] == "metric" and METHOD == "neighbors":
                method = f"""
    dict(
        input='Dropdown', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value="euclidean", clearable=False, options=['cityblock', 'cosine', 'euclidean', 'l1', 'l2', 'manhattan', 'braycurtis', 'canberra', 'chebyshev', 'correlation', 'dice', 'hamming', 'jaccard', 'kulsinski', 'mahalanobis', 'minkowski', 'rogerstanimoto', 'russellrao', 'seuclidean', 'sokalmichener', 'sokalsneath', 'sqeuclidean', 'yule'])
    )"""
            elif arg[0] == "var_names" and METHOD == "dendrogram":
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value={value},type="text")
    )"""
            elif arg[0] == "init_pos" and METHOD == "draw_graph":
                method = f"""
    dict(
        input='Dropdown', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value={value}, clearable=True, options=dict(function="['paga']+[i for i,j in config.adata.obsm.items() if j.shape[1] == 2]"))
    )"""
            elif arg[0] == "components" and METHOD == "embedding_density":
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value={value},type="text")
    )"""
            elif arg[0] == "key_added":
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description='{res}', 
        properties=dict(value='{value}',type="text")
    )"""
            elif arg[0] == "verbose":
                method = None
            elif arg[0] == "copy":
                method = None
            elif int == arg[2] or float == arg[2] or typing.Optional[int] == arg[2] or  typing.Optional[float] == arg[2]  or  typing.Union[float,int] == arg[2]:
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value={value},type="number")
    )"""
            elif bool == arg[2] or typing.Optional[bool] == arg[2]:
                method = f"""
    dict(
        input='BooleanSwitch', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(on={value})
    )"""
            elif str == arg[2] or typing.Optional[str] == arg[2]:
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value='{value}')
    )"""
            elif typing._LiteralGenericAlias == type(arg[2]):
                l = str(arg[2]).split("Literal")[-1]
                method = f"""
    dict(
        input='Dropdown', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value='{value}',options={l})
    )"""
            elif str(arg[2]).startswith("typing.Optional[typing.Literal["):
                l = str(arg[2]).split("Literal")[-1][:-1]
                method = f"""
    dict(
        input='Dropdown', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value='{value}',options={l})
    )"""
            elif int == type(arg[1]) or float == type(arg[1]):
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(value='{value}',type="number")
    )"""
            elif bool == type(arg[1]):
                method = f"""
    dict(
        input='BooleanSwitch', 
        name='{arg[0]}', 
        description="{res}", 
        properties=dict(on={value})
    )"""
            else:
                COUNT += 1
                method = None
                txt_ignored += "\n\t\t"+str(arg)
    #             method = f"""
    # dict(
    #     input='Input', 
    #     name='{arg[0]}', 
    #     description="{res}", 
    #     properties=dict(value='{value}')
    # )"""

            if method != None:
                executioncode += method+","
                codeargs += codearg

    if txt_missing != "" or txt_ignored != "":
        print(METHOD)
    if txt_missing != "":
        print("\tMissing:",txt_missing)
    if txt_ignored != "":
        print("\tIgnored:",txt_ignored)
    if txt_missing != "" or txt_ignored != "":
        print("\n")

    return executioncode, codeargs, n_args

s = dict()
for module,i in models:

    METHOD = i

    exec(f"l = [k for k in dir(sc.{module}.{i}) if k.startswith('get_')]",locals(),globals())
    for k in l:
        if k in s.keys():
            s[k].append(i)
        else:
            s[k] = [i]

    code = f"args = inspect.getfullargspec(sc.{module}.{i})"
    exec(code,locals(),globals())

    code_info = f"args_info = sc.{module}.{i}.__doc__"
    exec(code_info,locals(),globals())

    c, kargs, n_args = adapt(args, args_info)
    if n_args < 1:
        executioncode = "[ARGINPUT,"
        executioncode += c
        executioncode += "]"

        code = f"""
import numpy as np
from numpy import inf
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scanpy as sc

from general import *

{i.lower()}_args = dict(
    execution = {executioncode},
    postexecution = [],
    plot = []
)

def {i.lower()}_f(adata,kwargs):

    sc.{module}.{i}(
        adata,{kargs}
    )
        
    return

def {i.lower()}_plot():
    return []

config.methods["{i.lower()}"] = dict(
        
    properties=dict(
        type="QC",
        make_new_h5ad=False,
    ),

    args = {i.lower()}_args,

    function = {i.lower()}_f,

    plot = {i.lower()}_plot,

)
"""
        # exec(code, locals(), globals())
        file = f"../methods/{i.lower()}.py"
        with open(file,"w") as outfile:
            outfile.write(code)

print("Uncomplete functions: ", COUNT)