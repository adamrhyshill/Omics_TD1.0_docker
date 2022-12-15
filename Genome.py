from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import numpy as np
import pathlib
from app import app
from dash import dash_table
from Bio.KEGG.REST import kegg_find

PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../datasets").resolve()

countData = pd.read_csv(DATA_PATH.joinpath("R_Landing.csv")) 
countData = countData[['Fasta','Name','ECNumber','KEGGNAME','KEGGSYMBOL','GONAME_1','GONAME_2','AAs_sequence','Nucleotide_sequence']].copy()
countData = countData.rename(columns={'GONAME_1':'GONAME','GONAME_2':'GONAME 2','AAs_sequence':'AA','Nucleotide_sequence':'NS','ECNumber':'EC'})
countData = countData.sort_values('EC',ascending=False)
#countData['id'] = countData['Fasta']
#countData.set_index('id', inplace=True, drop=False)

layout = html.Div([
    dbc.Row(dbc.Col(html.H1("Genome Annotation Search", style={'textAlign':'center','fontSize': 30}),),),
    html.Hr(),
    dbc.Row(dbc.Col(dash_table.DataTable(
        id='halo-table-data',
        columns=[
            {"name": i, "id": i,
           #  "deletable": True,
              "selectable": True, "hideable": True}
            if i == "GONAME" or i == "GONAME 2" or i == "KEGGNAME" or i == "EC" or i == "AA" or i == "NS"
            else {"name": i, "id": i, 
            #"deletable": True, 
            "selectable": True}
            for i in countData.columns
        ],
        data=countData.to_dict('records'),  # the contents of the table
        editable=False,              # allow editing of data inside all cells
        filter_action="native",     # allow filtering of data by user ('native') or not ('none')
        sort_action="native",       # enables data to be sorted per-column by user or not ('none')
        sort_mode="multi",         # sort across 'multi' or 'single' columns
        column_selectable="multi",  # allow users to select 'multi' or 'single' columns
        row_selectable="multi",     # allow users to select 'multi' or 'single' rows
        row_deletable=False,         # choose if user can delete a row (True) or not (False)
        selected_columns=[],        # ids of columns that user selects
        selected_rows=[],           # indices of rows that user selects
        page_action="native",       # all data is passed to the table up-front or not ('none')
        page_current=0,             # page number that user is on
        page_size=10,                # number of rows visible per page
        style_cell={                # ensure adequate header width when text is shorter than cell's text
            'minWidth': 95, 'maxWidth': 95, 'width': 95
        },
        #style={'fontSize': 14},
        style_header={'backgroundColor': 'rgb(247, 211, 207)','color': 'black','fontWeight': 'bold','fontSize': 14, 'fontType':'Times'},
        style_cell_conditional=[    # align text columns to left. By default they are aligned to right
            {
                #'if': {'column_id': c},
                'textAlign': 'left', 'fontSize': 12, 'fontType':'Times'
            } #for c in ['country', 'iso_alpha3']
        ],
        style_data={  # overflow cells' content into multiple lines
            'whiteSpace': 'normal',
            'height': 'auto'
        },
    ) ,)),

    html.Br(),
    html.Br(),

    dbc.Row(dbc.Col(html.H3("Search the KEGG gene database, input 1-3 words for best query"))),
    dbc.Row(dbc.Col([dcc.Textarea(id='text_in',style={'width': '100%', 'height': 100, 'fontSize':14},),
    html.Button('Query', id='textarea-state-example-button', n_clicks=0)]),),
    dbc.Row(dbc.Col(dcc.Markdown(id='query_out', style={'width': '100%', 'height': 100, 'fontSize':14}))),
])


@app.callback(
    Output(component_id='query_out', component_property='children'),
    Input('textarea-state-example-button', 'n_clicks'),
    State('text_in', 'value'),
    )

def query_the_web(n_clicks,text_in):
    if n_clicks > 0:
        text_in = text_in.replace(" ",'+')
        q = kegg_find("genes", text_in).read()
        q = ' - ' + q.replace("\n", "\n - ")
        return q
