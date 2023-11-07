import inspect
import scanpy as sc
import typing
import re

models = [("pp",i) for i in dir(sc.pp) if not i.startswith("_") and i not in ["base","utils"]]
models += [("tl",i) for i in dir(sc.tl) if not i.startswith("_")  and i not in ["base","utils"]]
# models += [("epp",i) for i in dir(sc.external.pp) if not i.startswith("_") and i not in ["base","utils"]]
# models += [("etl",i) for i in dir(sc.external.tl) if not i.startswith("_") and i not in ["base","utils"]]

# for model in models:
#     print(model)

def adapt(args, args_info):

    try:
        largs = len(args.args)-len(args.defaults)
    except:
        largs = 0

    executioncode = ""
    kargs = ""
    for j,arg in enumerate(args.args):
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
            res = ""
            if j <= largs:
                value = None
            elif args.defaults == None:
                value = None
            else:
                value = args.defaults[j-largs]

            add = ""
            print(arg)
            if arg in args.annotations.keys():
                print(arg," ",args.annotations[arg])
                kargs += f"{arg}=kwargs['{arg}'],"
                if arg == "batch_key":
                    method = """
    dict(
        input = "Dropdown",
        name = "batch_key",
        description = "(Optional[str] (default: None)) - key in adata.obs for batch information. Categories will automatically be converted into integer categories and saved to adata.obs['_scanpy_batch']. If None, assigns the same batch to all the data.",
        properties = dict(
            value = None,
            clearable = True,
            options = dict(function = "[i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in ['str','int','categorical']]")
        )
    )"""
                elif arg == "layer":
                    method = """
    dict(
        input = "Dropdown",
        name = "layer",
        description = "(Optional[str] (default: None)) - if not None, uses this as the key in adata.layers for raw count data.",
        properties = dict(
            value = None,
            clearable = True,
            options = dict(function = "['X']+[i for i in config.adata.layers.keys()]")
        )
    )"""
                elif arg == "verbose":
                    None
                elif arg == "copy":
                    None
                elif int == args.annotations[arg] or float == args.annotations[arg]:
                    method = f"""
    dict(
        input='Input', 
        name='{arg}', 
        description="{res}", 
        properties=dict(value={value},type="number")
    )"""
                elif bool == args.annotations[arg]:
                    method = f"""
    dict(
        input='BooleanSwitch', 
        name='{arg}', 
        description="{res}", 
        properties=dict(on={value})
    )"""
                elif str == args.annotations[arg]:
                    method = f"""
    dict(
        input='Input', 
        name='{arg}', 
        description="{res}", 
        properties=dict(value='{value}')
    )"""
                elif typing._LiteralGenericAlias == type(args.annotations[arg]):
                    l = str(args.annotations[arg]).split("Literal")[-1]
                    method = f"""
    dict(
        input='Dropdown', 
        name='{arg}', 
        description="{res}", 
        properties=dict(value='{value}',options={l})
    )"""
                else:
                    method = f"""
    dict(
        input='Input', 
        name='{arg}', 
        description="{res}", 
        properties=dict(value='{value}')
    )"""

                add += method+""
                executioncode += add+","

    return executioncode, kargs

s = dict()
for module,i in models:

    exec(f"l = [k for k in dir(sc.{module}.{i}) if k.startswith('get_')]",locals(),globals())
    for k in l:
        if k in s.keys():
            s[k].append(i)
        else:
            s[k] = [i]

    code = f"args = inspect.getfullargspec(sc.{module}.{i})"
    exec(code,locals(),globals())

    print(args)

    code_info = f"args_info = sc.{module}.{i}.__doc__"
    exec(code_info,locals(),globals())

    executioncode = "[ARGINPUT,"
    c, kargs = adapt(args, args_info)
    executioncode += c
    executioncode += "]"

    code = f"""
import numpy as np
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

    sc.{module}.{i}(adata,{kargs})
    
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

    file = f"scanpy_methods/{i.lower()}.py"
    with open(file,"w") as outfile:
        outfile.write(code)