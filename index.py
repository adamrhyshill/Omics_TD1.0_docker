import dash
from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output

# Connect to main app.py file
# imports any asset objects into the page. if you add anything to app.py, add here
from app import app
from app import server

# Connect to your app pages; import all the python code
from apps import KEGG_map, PCA, Genome#, Flux

sidebar = dbc.Card([dbc.CardBody(
    [
        html.H2("Halomonas TD1.0 OMICs", className="display-5"),
        html.Hr(),
        html.P(
            "Analysis Options", className="lead",style={'fontSize': 14}
        ),
        dbc.Nav(
            [
                dbc.NavLink("PCA", href='/apps/PCA', active="exact",style={'fontSize': 14}),
                dbc.NavLink("KEGG Map", href='/apps/KEGG_map', active="exact",style={'fontSize': 14}),
              #  dbc.NavLink("GEM Flux", href='/apps/Flux', active="exact",style={'fontSize': 14}),
                dbc.NavLink("Genome", href='/apps/Genome', active="exact",style={'fontSize': 14}),
            ],
            vertical=True,
            pills=True,
        ),
    ]),], color='light', style={"height":"100vh", "width":"18rem","position":"fixed"})


content = html.Div(id="page-content", children=[], style={"padding": "0.5rem"},)

app.layout = dbc.Container([
    dcc.Location(id="url", refresh=False),
    dbc.Row([dbc.Col(sidebar,width=2), dbc.Col(content,width=10,style={"margin-left":"18rem"})
    ])
],fluid=True)

@app.callback(
    Output("page-content", "children"),
    [Input("url", "pathname")])

def render_page_content(pathname, value = '/apps/PCA'):
    if pathname == '/apps/KEGG_map':
        return KEGG_map.layout
    elif pathname == '/apps/PCA':
        return PCA.layout
  #  elif pathname == '/apps/Flux':
  #      return Flux.layout
    elif pathname == '/apps/Genome':
        return Genome.layout
    else:
        return html.Div(
        [
            html.H1("Welcome to Halomonas TD1.0 OMICS landing page", style={'textAlign':'center','fontSize': 30}),
            html.Hr(),
            dcc.Markdown(f"The pathname **{pathname}** was not recognised. Please select an analysis option at left.", style={'textAlign':'center','fontSize': 12, 'color':'red'}),
        ]
    )

if __name__ == '__main__':
    app.run_server(
        host='127.0.0.1',port='8043',
        debug=True, )
