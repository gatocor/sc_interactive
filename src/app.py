import dash
import dash_bootstrap_components as dbc
import dash_cytoscape as cyto
cyto.load_extra_layouts()

# bootstrap theme
# https://bootswatch.com/lux/
external_stylesheets = [dbc.themes.LUX]

app = dash.Dash(__name__, external_stylesheets=external_stylesheets) # Local stylesheet (LUX)
# app = dash.Dash(__name__, external_stylesheets=external_stylesheets) # External stylesheet

server = app.server
app.config.suppress_callback_exceptions = True