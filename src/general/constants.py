ENDH5AD = ".h5ad"
ENDANALYSIS = ".sc"

RAWNODE = {
        'data': {
            'id':'Raw', 
            'name':'Raw',
            'type':'Raw', 
            'method':'Raw', 
            'label':'Raw', 
            'color':'white',
            'h5ad_file':'Raw.h5ad',
            'computed':True,
            'opacity':1,
            'summary':'Raw', 
            'image':'../assets/Raw.png',
            'parameters':{'input':None}}, 
        'position':{'x':0,'y':0},
    }

STARTPATH = '../../'  # Replace with the actual path to the folder containing .h5ad files

ARGINPUT = {
            "input":"Input",
            "name":"input",
            "description":"Incoming analysis node.",
            "properties":{
                "value": None,
                "readonly": True
            }
        }

ARGMATRIX = {
            "input":"Dropdown",
            "name":"matrix",
            "description":"Matrix from the adata object to use (to choose between X, layers or obsm; if the two later have any key).",
            "properties":{
                "value":"X",
                "clearable":False,
                "options": {"function":"matrix_options()"}
            }
        }

ARGMATRIX_KEY = {
            "input":"Dropdown",
            "name":"matrix_key",
            "description":"Observable matrix to use",
            "properties":{
                "value": {"function":"matrix_keys()[0]"},
                "clearable":False,
                "options": {"function":"matrix_keys()"},
            },
            "visible":{"function":"config.active_node_parameters['matrix'] in ['layers','obsm']"}
        }

ARGBATCH = {
            "input":"Dropdown",
            "name":"batch",
            "description":"adada.obs column used to perform the current node analysis in batches by this key.",
            "properties":{
                "value":None,
                "clearable":True,
                "options": {"function":"[i for i,j in zip(config.adata.obs.columns.values, config.adata.obs.dtypes) if j in ['str','object','category','int']]"}
            }
        }

ARGS_REPRESENTATION = [
        {
            "input":"Dropdown",
            "name":"plot_representation",
            "description":"Representation used to plot.",
            "properties":{
                "value":{"function":"[i for i in config.adata.obsm.keys()][0]"},
                "clearable":True,
                "options":{"function":"[i for i in config.adata.obsm.keys()]"} 
            },
        },
        {
            "input":"Dropdown",
            "name":"plot_style",
            "description":"Representation used to plot.",
            "properties":{
                "value":"scatter_2d",
                "clearable":False,
                "options":{"function":"get_representation_styles()"} 
            },
            "recomputeAfter":["plot_representation"]
        },
        {
            "input":"Dropdown",
            "name":"plot_n_components",
            "description":"Number of component to plot.",
            "properties":{
                "value":2,
                "clearable":True,
                "options":{"function":"np.arange(2,config.adata.obsm[get_node(config.selected)['data']['plot']['plot_representation']].shape[1]+1)"} 
            },
            "visible":"get_node(config.selected)['data']['plot']['plot_style'] == 'scatter_2d'",
            "recomputeAfter":["plot_representation"]
        },
        {
            "input":"Dropdown",
            "name":"plot_dimension_x",
            "description":"Dimension to plot in x axis.",
            "properties":{
                "value":0,
                "clearable":True,
                "options":{"function":"np.arange(0,config.adata.obsm[get_node(config.selected)['data']['plot']['plot_representation']].shape[1])"} 
            },
            "visible":"get_node(config.selected)['data']['plot']['plot_style'] in ['scatter_2d','scatter_3d']",
            "recomputeAfter":["plot_representation"]
        },
        {
            "input":"Dropdown",
            "name":"plot_dimension_y",
            "description":"Dimension to plot in x axis.",
            "properties":{
                "value":1,
                "clearable":True,
                "options":{"function":"np.arange(0,config.adata.obsm[get_node(config.selected)['data']['plot']['plot_representation']].shape[1])"} 
            },
            "visible":"get_node(config.selected)['data']['plot']['plot_style'] in ['scatter_2d','scatter_3d']",
            "recomputeAfter":["plot_representation"]
        },
        {
            "input":"Dropdown",
            "name":"plot_dimension_z",
            "description":"Dimension to plot in x axis.",
            "properties":{
                "value":2,
                "clearable":True,
                "options":{"function":"np.arange(0,config.adata.obsm[get_node(config.selected)['data']['plot']['plot_representation']].shape[1])"} 
            },
            "visible":"get_node(config.selected)['data']['plot']['plot_style'] in ['scatter_3d']",
            "recomputeAfter":["plot_representation"]
        },
]

ARGS_COLOR = [
        {
            "input":"Dropdown",
            "name":"color",
            "description":"Color used to represent.",
            "properties":{
                "value":None,
                "clearable":True,
                "options":{"function":"args_color()"} 
            },
        },
        {
            "input":"Dropdown",
            "name":"color_var",
            "description":"What var key to use.",
            "properties":{
                "value":{"function":"args_color_var()[0]"},
                "clearable":False,
                "options":{"function":"args_color_var()"} 
            },
            "visible":{"function":"None != args_color_var()[0]"},
            "recomputeAfter": ["color"] 
        },
        {
            "input":"Dropdown",
            "name":"color_layer",
            "description":"What var key to use.",
            "properties":{
                "value":"X",
                "clearable":False,
                "options":{"function":"['X']+[i for i in config.adata.layers.keys()]"} 
            },
            "visible":{"function":"None != args_color_var()[0]"},
            "recomputeAfter": ["color"] 
        },
        {
            "input":"Dropdown",
            "name":"color_obsm_dimension",
            "description":"What var key to use.",
            "properties":{
                "value":{"function":"args_color_obsm()[0]"},
                "clearable":False,
                "options":{"function":"args_color_obsm()"} 
            },
            "visible":{"function":"None != args_color_obsm()[0]"},
            "recomputeAfter": ["color"] 
        },
]

GRAPHSTYLESHEET = [
                    {
                        'selector': 'node',
                        'style': {
                            'label': 'data(id)',  # Show the node labels
                            'text-wrap': 'wrap',
                            # 'text-max-width': 30,
                            # 'text-justification': 'left',
                            'text-margin-x':-21,
                            'text-margin-y':-24,
                            'text-valign': 'bottom',  # Center the label vertically
                            'text-halign': 'right',  # Center the label horizontally
                            'background-color': 'white',
                            'opacity':'data(opacity)',
                            # 'shape': 'box',
                            'width': 24.5,
                            'height': 25,
                            'border-color': 'black',
                            'border-width': .0,
                            'shape': 'square',
                            'font-size':2,
                            'z-index-compare':'manual',
                            'z-index':0,
                            'background-fit': 'cover',
                            'background-image': 'data(image)'
                        },
                    },
                    {
                        'selector': 'edge',
                        'style': {
                            'curve-style':'straight',
                            'line-cap':'round',
                            'source-endpoint':'outside-to-node',
                            'source-arrow-shape': 'circle',
                            'target-arrow-shape': 'circle',
                            'mid-target-arrow-shape': 'triangle',
                            'width':1,
                            'arrow-scale': .35,
                            'z-index-compare':'manual',
                            'z-index':1,
                            'source-endpoint': '10.5 0', 
                            'target-endpoint': '-10.5 0',
                            # 'line-color':'black'
                        },
                    },
                    {
                        'selector': ':selected',
                        'style': {
                            'border-color': 'red'
                        }
                    },
                ]

GRAPHLAYOUT = {'name':'dagre',
                'roots':'[id="Raw"]',
                'rankDir': "LR",
                'rankSep': 5}

ARGUMENTBOXSTYLE = {"background-color":"lightgray"}