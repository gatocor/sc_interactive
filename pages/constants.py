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

STARTPATH = '../'  # Replace with the actual path to the folder containing .h5ad files

ARGINPUT = {
            "input":"Dropdown",
            "name":"input",
            "description":"Observable to use",
            "value":None,
            "clearable":False,
            "options": {"function": "node_names(exclude_downstream_from_node=config.selected)"}
        }

ARGBATCH = {
            "input":"Dropdown",
            "name":"batch",
            "description":"Observable to use",
            "value":None,
            "clearable":True,
            "options": {"function":"[i for i,j in zip(config.adata.obs.columns.values, config.adata.obs.dtypes) if j in ['str','object','category','int']]"}
        }

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

ARGUMENTBOXSTYLE = {"background-color":"lightgray"}