import scvi
import typing
import re

models = [("model",i) for i in dir(scvi.model) if not i.startswith("_") and i not in ["base","utils"]]
models += [("external",i) for i in dir(scvi.external) if not i.startswith("_") and i[0].isupper() and i not in ["base","utils"]]

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
            if arg in args.annotations.keys():
                kargs += f"{arg}=kwargs['{arg}'],"
                if arg == "batch_key":
                    method = """
    dict(
        input = "Dropdown",
        name = "batch_key",
        description = "(Optional[str] (default: None)) - key in adata.obs for batch information. Categories will automatically be converted into integer categories and saved to adata.obs['_scvi_batch']. If None, assigns the same batch to all the data.",
        properties = dict(
            value = None,
            clearable = True,
            options = dict(function = "[i for i in config.adata.obs.columns.values if config.adata.obs.dtypes[i] in ['str','int','categorical']]")
        )
    )"""
                elif arg == "labels_key":
                    method = """
    dict(
        input = "Dropdown",
        name = "labels_key",
        description = "(Optional[str] (default: None)) - key in adata.obs for label information. Categories will automatically be converted into integer categories and saved to adata.obs['_scvi_labels']. If None, assigns the same label to all the data.",
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
                elif arg == "n_hidden":
                    method = """
    dict(
        input = "Dropdown",
        name = "n_hidden",
        description = "(int (default: 128)) - Number of nodes per hidden layer.",
        properties = dict(
            value = 128,
            clearable = False,
            options = [2**i for i in range(10)]
        )
    )"""
                elif arg == "n_latent":
                    method = """
    dict(
        input = "Dropdown",
        name = "n_latent",
        description = "(int (default: 10)) - Dimensionality of the latent space.",
        properties = dict(
            value = 10,
            clearable = False,
            options = [i for i in range(100)]
        )
    )"""
                elif arg == "n_layers":
                    method = """
    dict(
        input = "Dropdown",
        name = "n_layers",
        description = "(int (default: 1)) - Number of hidden layers used for encoder and decoder NNs.",
        properties = dict(
            value = 1,
            clearable = False,
            options = [i for i in range(100)]
        )
    )"""
                elif arg == "dropout_rate":
                    method = """
    dict(
        input = "Input",
        name = "dropout_rate",
        description = "(float (default: 0.1)) - Dropout rate for neural networks.",
        properties = dict(
            value = .1,
            type = "number"
        )
    )"""
                elif arg == "max_epochs":
                    method = """
    dict(
        input = "Input",
        name = "max_epochs",
        description = "(Optional[int] (default: None)) - Number of passes through the dataset. If None, defaults to np.min([round((20000 / n_cells) * 400), 400])",
        properties = dict(
            value = None,
            type = "number"
        ),
    )"""
                elif arg == "accelerator":
                    method = """
    dict(
        input = "Dropdown",
        name = "accelerator",
        description = "(str (default: 'auto')) - Supports passing different accelerator types ('cpu', 'gpu', 'tpu', 'ipu', 'hpu', 'mp', 'auto') as well as custom accelerator instances.",
        properties = dict(
            value = "auto",
            options = ['cpu', 'gpu', 'tpu', 'ipu', 'hpu', 'mp', 'auto']
        ),
    )"""
                elif arg == "devices":
                    method = """
    dict(
        input = "Input",
        name = "devices",
        description = "(Union[int, List[int], str] (default: 'auto')) - The devices to use. Can be set to a non-negative index (int or str), a sequence of device indices (list or comma-separated str), the value -1 to indicate all available devices, or 'auto' for automatic selection based on the chosen accelerator. If set to 'auto' and accelerator is not determined to be 'cpu', then devices will be set to the first available device.",
        properties = dict(
            value = "auto"
        ),
    )"""
                elif arg == "train_size":
                    method = """
    dict(
        input = "Input",
        name = "train_size",
        description = "(float (default: 0.9)) - Size of training set in the range [0.0, 1.0].",
        properties = dict(
            value = 0.9,
            type = "number"
        ),
    )"""
                elif arg == "validation_size":
                    method = """
    dict(
        input = "Input",
        name = "validation_size",
        description = "(Optional[float] (default: None)) - Size of the test set. If None, defaults to 1 - train_size. If train_size + validation_size < 1, the remaining cells belong to a test set.",
        properties = dict(
            value = .1,
            type = "number"
        ),
    )"""
                elif arg == "shuffle_set_split":
                    method = """
    dict(
        input = "BooleanSwitch",
        name = "shuffle_set_split",
        description = "(bool (default: True)) - Whether to shuffle indices before splitting. If False, the val, train, and test set are split in the sequential order of the data according to validation_size and train_size percentages.",
        properties = dict(
            on = True
        ),
    )"""
                elif arg == "batch_size":
                    method = """
    dict(
        input = "Dropdown",
        name = "batch_size",
        description = "(int (default: 128)) - Minibatch size to use during training.",
        properties = dict(
            value = 128,
            options = [2**i for i in range(10)]
        ),
    )"""
                elif arg == "early_stopping":
                    method = """
    dict(
        input = "BooleanSwitch",
        name = "early_stopping",
        description = "(bool (default: False)) - Perform early stopping. Additional arguments can be passed in **kwargs. See Trainer for further options.",
        properties = dict(
            on = False
        ),
    )"""
                elif arg == "plan_kwargs":
                    None
                elif int == args.annotations[arg] or float == args.annotations[arg]:
                    method = f"""
    dict(
        input='Input', 
        name='{arg}', 
        description="{res}", 
        properties=dict(value={value},type="number")
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

    exec(f"l = [k for k in dir(scvi.{module}.{i}) if k.startswith('get_')]",locals(),globals())
    for k in l:
        if k in s.keys():
            s[k].append(i)
        else:
            s[k] = [i]

    code_prepare = f"args_prepare = inspect.getfullargspec(scvi.{module}.{i}.setup_anndata)"
    exec(code_prepare,locals(),globals())

    code_prepare_info = f"args_prepare_info = scvi.{module}.{i}.setup_anndata.__doc__"
    exec(code_prepare_info,locals(),globals())

    code = f"args = inspect.getfullargspec(scvi.{module}.{i})"
    exec(code,locals(),globals())

    code_info = f"args_info = scvi.{module}.{i}.__doc__"
    exec(code_info,locals(),globals())

    code_train = f"args_train = inspect.getfullargspec(scvi.{module}.{i}.train)"
    exec(code_train,locals(),globals())

    code_train_info = f"args_train_info = scvi.{module}.{i}.train.__doc__"
    exec(code_train_info,locals(),globals())

    executioncode = "[ARGINPUT,"
    c, kargs =  adapt(args_prepare, args_prepare_info)
    executioncode += c
    c, kargs2 = adapt(args, args_info)
    executioncode += c
    c, kargs3 = adapt(args_train, args_train_info)
    executioncode += c
    executioncode += "]"

    code = f"""
import numpy as np
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import scvi
scvi.settings.seed = 0

from general import *

{i.lower()}_args = dict(
    execution = {executioncode},
    postexecution = [],
    plot = []
)

def {i.lower()}_f(adata,kwargs):

    scvi.model.{i}.setup_anndata(adata,{kargs})
    model = scvi.{module}.{i}(adata,{kargs2})
    model.train({kargs3})
    
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

    file = f"scvi_methods/{i.lower()}.py"
    with open(file,"w") as outfile:
        outfile.write(code)

for i,j in s.items():
    print(i,j)