import inspect
import scanpy as sc
import typing
import re

include = [
    'calculate_qc_metrics', 
    'combat', 
    'downsample_counts', 
    'filter_cells', 
    'filter_genes', 
    'filter_genes_dispersion', 
    'highly_variable_genes', 
    'log1p', 
    'neighbors', 
    'normalize_total', 
    'pca', 
    'recipe_seurat', 
    'recipe_weinreb17', 
    'recipe_zheng17', 
    'regress_out', 
    'scale', 
    'sqrt', 
    'subsample', 
    'dendrogram', 
    'diffmap', 
    'dpt', 
    'draw_graph', 
    'embedding_density', 
    'filter_rank_genes_groups', 
    'ingest', 
    'leiden', 
    'louvain', 
    'marker_gene_overlap', 
    'paga', 
    'paga_compare_paths', 
    'paga_degrees', 
    'paga_expression_entropies', 
    'pca', 
    'rank_genes_groups', 
    'score_genes', 
    'score_genes_cell_cycle', 
    'tsne', 
    'umap'         
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
POSITION = "parameters"

def adapt(args, args_info):

    global COUNT, METHOD, POSITION

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
        try:
            res = re.findall(arg+r'\n.*\n', args_info)[0].split("\n")[1].split("  ")[-1]
        except:
            res = ""
        res = res.replace('`','\'')
        res = res.replace('\'"','\'')
        res = res.replace('"\'','\'')
        res = arg[2]

        t = re.findall(r'<class \'.*?\'>', str(arg[2]))
        if len(t) > 0:
            t = t[0].split('\'')[1]
        else:
            t = arg[2]

        if arg[2] != None:
            if isinstance(arg[1],str):
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        visible=dict(function="'{arg[1]}'!=eval(get_node(config.selected)['data']['{POSITION}']['{arg[0]}']) or config.show_{POSITION}"),
        properties=dict(value="'{arg[1]}'",type="text")
    )"""
                codearg = f"""
        {arg[0]}=type_formater(kwargs["{arg[0]}"],{t}),"""

            else:
                method = f"""
    dict(
        input='Input', 
        name='{arg[0]}', 
        description="{res}", 
        visible=dict(function="str({arg[1]})!=get_node(config.selected)['data']['{POSITION}']['{arg[0]}'] or config.show_{POSITION}"),
        properties=dict(value="{arg[1]}",type="text")
    )"""
                codearg = f"""
        {arg[0]}=type_formater(kwargs["{arg[0]}"],{t}),"""

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
    POSITION = "parameters"

    exec(f"l = [k for k in dir(sc.{module}.{i}) if k.startswith('get_')]",locals(),globals())
    for k in l:
        if k in s.keys():
            s[k].append(i)
        else:
            s[k] = [i]

    #Execution
    code = f"args = inspect.getfullargspec(sc.{module}.{i})"
    exec(code,locals(),globals())
    code_info = f"args_info = sc.{module}.{i}.__doc__"
    exec(code_info,locals(),globals())
    c, kargs, n_args = adapt(args, args_info)

    ps = [j for j in plot if i in j[1]]
    figs = ""
    plot_c = ""
    POSITION = "plot"
    if len(ps)>0:
        for plot_module, plot_function in ps[:1]:
            code = f"args = inspect.getfullargspec(sc.{plot_module}.{plot_function})"
            exec(code,locals(),globals())
            code_info = f"args_info = sc.{plot_module}.{plot_function}.__doc__"
            exec(code_info,locals(),globals())
            plot_c, plot_kargs, plot_n_args = adapt(args, args_info)

            figs += f"""
    fig = sc.{plot_module}.{plot_function}(
        config.adata,{plot_kargs}
    )
"""

    if n_args <= 1:
        executioncode = "[ARGINPUT,"
        executioncode += c
        executioncode += "]"

        plot_executioncode = "["
        plot_executioncode += plot_c
        plot_executioncode += "]"


        code = f"""
import numpy
from numpy import inf
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scanpy as sc
import louvain
import scipy
import leidenalg
import plotly.tools as tls
import cycler
import matplotlib      # pip install matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import base64
from io import BytesIO
from general import *

{i.lower()}_args = dict(
    execution = {executioncode},
    postexecution = [],
    plot = {plot_executioncode}
)

def {i.lower()}_f(adata,kwargs):

    sc.{module}.{i}(
        adata,{kargs}
    )
        
    return

def {i.lower()}_plot():

    kwargs = get_node(config.selected)['data']['plot']
    {figs}

    # Save it to a temporary buffer.
    buf = BytesIO()
    fig.savefig(buf, format="png")
    # Embed the result in the html output.
    fig_data = base64.b64encode(buf.getbuffer()).decode("ascii")
    fig_bar_matplotlib = f'data:image/png;base64,'+fig_data
    fig =  html.Img(id='bar-graph-matplotlib',src=fig_bar_matplotlib)

    return fig

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
        
        code = code.replace("NoneType","str")
        code = re.sub(r'return_fig=.*,',"return_fig=True,", code)
        # exec(code, locals(), globals())
        file = f"../methods/{i.lower()}.py"
        with open(file,"w") as outfile:
            outfile.write(code)

print("Uncomplete functions: ", COUNT)