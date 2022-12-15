from dash import dcc
#import dash
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Input, Output
import plotly.express as px
#import plotly.graph_objects as go
import pandas as pd
import numpy as np
from sklearn import decomposition
#from dash import dash_table
from dash_bio import Clustergram as clust
import pathlib
from app import app


# -------------------------------------------------------------------------------------
# functions

def PCA_script(df_genes, transpose = False, loadings = False, pcs=2): #PCA with 2 pc's
    matrix = df_genes.to_numpy()
    if transpose:
        matrix = matrix.T
    pca = decomposition.PCA(n_components=pcs) #run decomposision
    pca.fit(matrix)
    if loadings: #returns to analyze and also returns the loadings for each variable
        return(pca.transform(matrix), pca.components_.T * np.sqrt(pca.explained_variance_), pca.explained_variance_ratio_)
    return(pca.transform(matrix))

def standardize_2(df):
    df = df[df>0]
    df = np.log2(df)
    df = (df - df.mean())/df.std()
    return df

def make_df(choice, PCA):

    if choice == 'all_data': 
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLNaCL_ProteinpgDW', '20gLNaCL_ProteinpgDW','100gLNaCL_ProteinpgDW','Urea3.6gL_ProteinpgDW','09Fermentation_ProteinpgDW', '19Fermentation_ProteinpgDW','30Fermentation_ProteinpgDW', '20gLNaCL_RPKM', '60gLNaCL_RPKM','100gLNaCL_RPKM','Urea3.6gL_RPKM','09Fermentation_RPKM','19Fermentation_RPKM', '30Fermentation_RPKM']].copy()
    elif choice == 'fermentation_protein':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','09Fermentation_ProteinpgDW', '19Fermentation_ProteinpgDW', '30Fermentation_ProteinpgDW']].copy()
    elif choice == 'fermentation_transcript':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','09Fermentation_RPKM','19Fermentation_RPKM', '30Fermentation_RPKM']].copy()
    elif choice == 'fermentation_all':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','09Fermentation_RPKM','19Fermentation_RPKM', '30Fermentation_RPKM','09Fermentation_ProteinpgDW', '19Fermentation_ProteinpgDW','30Fermentation_ProteinpgDW']].copy()
    elif choice == 'shake_flask_protein':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLNaCL_ProteinpgDW', '20gLNaCL_ProteinpgDW','100gLNaCL_ProteinpgDW', 'Urea3.6gL_ProteinpgDW']].copy()
    elif choice == 'shake_flask_transcript':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','20gLNaCL_RPKM', '60gLNaCL_RPKM','100gLNaCL_RPKM','Urea3.6gL_RPKM']].copy()
    elif choice == 'all_transcript':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLNaCL_RPKM', '20gLNaCL_RPKM','100gLNaCL_RPKM','Urea3.6gL_RPKM','09Fermentation_RPKM','19Fermentation_RPKM', '30Fermentation_RPKM']].copy()
    elif choice == 'all_protein':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLNaCL_ProteinpgDW','20gLNaCL_ProteinpgDW','100gLNaCL_ProteinpgDW','Urea3.6gL_ProteinpgDW', '09Fermentation_ProteinpgDW', '19Fermentation_ProteinpgDW','30Fermentation_ProteinpgDW']].copy()        
    elif choice == 'shake_flask_all':
        PCA = PCA[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLNaCL_RPKM', '20gLNaCL_RPKM','100gLNaCL_RPKM', 'Urea3.6gL_RPKM','20gLNaCL_ProteinpgDW', '60gLNaCL_ProteinpgDW','100gLNaCL_ProteinpgDW', 'Urea3.6gL_ProteinpgDW']].copy()
    for column in PCA.columns.drop(['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL']):
        PCA[column] = standardize_2(PCA[column])
    PCA = PCA.dropna()
    return PCA

# -------------------------------------------------------------------------------------
# read in data # import data and format it
PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../datasets").resolve()

pathlist_2 = ['Ribosom', 'Carbon fixation pathways in prokaryotes', 'ABC transporter', 'Glyoxylate and dicarboxylate metabolism', 'Arginine biosynthesis', 'Butanoate metabolism', 'Quorum sensing',  'Fatty acid metabolism',  'Citrate cycle (TCA cycle)', 'Pyrimidine metabolism', '2-Oxocarboxylic acid metabolism', 'Synthesis and degradation of ketone bodies', 'Biosynthesis of amino acids','Glycine, serine and threonine metabolism', 'Protein export', 'Pentose phosphate pathway', 'Arginine and proline metabolism',  'Pyruvate metabolism', 'Oxidative phosphorylation','Phosphotransferase system (PTS)',  'DNA replication',   'Glutathione metabolism',    'Tyrosine metabolism', 'Lysine degradation','Homologous recombination', 'Folate biosynthesis', 'Galactose metabolism','Phenylalanine metabolism', 'RNA polymerase', 'Biofilm formation - Escherichia coli', 'Flagellar assembly',  'RNA degradation', 'C5-Branched dibasic acid metabolism', 'Bacterial chemotaxis', 'RNA transport', 'Vitamin B6 metabolism','Purine metabolism', 'Propanoate metabolism', 'Lipoic acid metabolism', 'Glycolysis / Gluconeogenesis','Selenocompound metabolism', 'Aminoacyl-tRNA biosynthesis','Amino sugar and nucleotide sugar metabolism','Lipopolysaccharide biosynthesis', 'Base excision repair','Nicotinate and nicotinamide metabolism', 'Ubiquinone and other terpenoid-quinone biosynthesis','Alanine, aspartate and glutamate metabolism', 'Glycerophospholipid metabolism', 'Thiamine metabolism',  'Peptidoglycan biosynthesis', 'Mismatch repair','Phenylalanine, tyrosine and tryptophan biosynthesis', 'Sulfur metabolism', 'Cysteine and methionine metabolism', 'Ascorbate and aldarate metabolism', 'Methane metabolism', 'Biosynthesis of unsaturated fatty acids','Terpenoid backbone biosynthesis', 'Riboflavin metabolism', 'Sulfur relay system',  'D-Glutamine and D-glutamate metabolism','Nucleotide excision repair', 'Thermogenesis', 'Pyrimidine metabolism', 'One carbon pool by folate', 'Nitrogen metabolism', 'Glycerolipid metabolism', 'Inositol phosphate metabolism','Glucosinolate biosynthesis',       'Fructose and mannose metabolism', 'Valine, leucine and isoleucine biosynthesis',  'Valine, leucine and isoleucine degradation', 'Fatty acid degradation']
countData = pd.read_csv(DATA_PATH.joinpath("R_Landing.csv"))  
PCA_data_options = ['all_data','all_transcript','all_protein','fermentation_all','shake_flask_all','shake_flask_protein','fermentation_protein','shake_flask_transcript','fermentation_transcript']
countData['KEGG_Pathway'] = 'Missing'
i=0
for key in countData.KEGGPATH:
    for kegginput in pathlist_2:
        if kegginput in key:
            countData.loc[i,'KEGG_Pathway'] = kegginput
    i=i+1

# -------------------------------------------------------------------------------------
# App layout

# define rows and THEN columns
# width defines how many columns the page gets. 6 columns gets 1/2 of the page!
    # if don't define width, it will expand to entire page

    # cheatsheet! https://dashcheatsheet.pythonanywhere.com

    # size (1-12, true (each element 1/2 page), auto (to size needs)),
    # offset: (1-12) number of invisible columns to left of component

    #'color': 'blue', for titles

    # CardGroup, CardDeck (both arrange in row nicely), and CardColumn(this one like pinterest)

    # colors https://plotly.com/python/discrete-color/
    # hex generator https://colordesigner.io/gradient-generator

card_content_right = html.Div([#dbc.CardHeader("Compounds",style={'fontSize': 16}),
    dbc.Button(
            "Info", style={'fontSize': 14},
            id="hover-target", 
           # color="danger",
            className="me-1",
            n_clicks=0,
        ),
        dbc.Popover(
            dcc.Markdown("Filter the cluster heatmap by **clicking** on a KEGG pathway in the bar chart. Refresh page to remove filter.", style={'fontSize': 14}), 
            target="hover-target",placement="bottom",offset="50,20",
            body=True,
            trigger="hover", 
        ),
    
     ])
card_content_left = [#dbc.CardHeader("Genes", style={'fontSize': 16}),
    dbc.Row([dbc.Col(html.H4("Select Principle Component, and pathway weight method (mean / sum).", style={'fontSize': 16}),
                        width={'size': 12, 'offset': 0}),
                  
                        
        ],className="mx-2"),
        dbc.Row([dbc.Col(dcc.Dropdown(id='loadings_for_graph', value='loadings PC2', clearable=False,
                                    persistence=True, persistence_type='memory',
                                    options=['loadings PC1','loadings PC2','loadings PC3','loadings PC4'],
                                    style={'fontSize': 14}),width={'size': 2, "offset": 0}),
            dbc.Col(dbc.RadioItems(id='av_or_sum',value='Average',options=[{'label':'Average','value':'Average'},{'label':'Sum','value':'Sum'}],
                                    style={'fontSize': 14}),width={'size': 1, "offset": 0}),
            dbc.Col(card_content_right,width={'size': 1, "offset": 0})
            
            
            ],className="mx-2" )]
info_content_right = html.Div([
    dbc.Button(
            "Info", style={'fontSize': 14},
            id="hover-target-2", 
           # color="danger",
            className="me-1",
            n_clicks=0,
        ),
        dbc.Popover(
            dcc.Markdown('''
            - PCA is an unsupervised machine-learning method that aims to reduce dimensionality. 
            
            - PCA can take in large data sets, and through linear transformation, simplify while minimizing information loss. 
            
            - Therefore, we can learn which indicators drive the most variance between different conditions.
            ''', 
            style={'fontSize': 12}), 
            target="hover-target-2",placement="bottom",offset="50,20",
            body=True,
            trigger="hover", 
        ),
    
     ])
    # dbc.Row([dbc.Col(dcc.Markdown("Filter the cluster heatmap by **clicking** on a KEGG pathway in the bar chart. Refresh page to remove filter.", style={'fontSize': 16}),
    #                     width={'size': 12, 'offset': 0}),
    #               ],className="card-text mx-2 pt-1"),
    

    #st.markdown("PCA is an unsupervised machine-learning method that aims to reduce dimensionality. PCA can take in large data sets, 
    # and through linear transformation, simplify while minimizing information loss. Therefore, we can learn which indicators drive the most variance between different conditions. [1]")


layout = html.Div([
        dbc.Row(dbc.Col(html.H1("PCA, Pathway Enrichments Analysis, and Heatmap", style={'textAlign':'center','fontSize': 30}),
                       
                        #width={'size': 6, 'offset': 3},md={'size':'true'},
                        ),),
        html.Hr(),
        html.Div(dbc.Row(dbc.Col(html.H4("Select OMICs and condition data:", style={'fontSize': 16}),
                        width={'offset': 0})),className="mx-2"),
        dbc.Row([dbc.Col(dcc.Dropdown(id='PCA_data_options', value='all_data', clearable=False,
                                    persistence=True, persistence_type='memory',
                                    options=[{'label': x, 'value': x} for x in PCA_data_options],
                                    style={'fontSize': 14}), 
                                    width={'size': 2, "offset":0}),
            dbc.Col(dbc.RadioItems(id='KEGG_missing_included',value='Include all data points',
                                    options=[
                                        {'label':'Include all data points','value':'Include all data points'},
                                        {'label':'Remove "Missing" KEGG Pathways','value':'Remove "Missing" KEGG Pathways'}],
                                    style={'fontSize': 14}),
                                    width={'size': 3, "offset": 0}),
            dbc.Col(info_content_right,width={'size': 1, "offset": 0}),
            ], className="mx-2"),
        html.Hr(),
        dbc.CardGroup(
            [
                dcc.Graph(id='PCA-chart', figure={},config={'showTips':True}),
                dcc.Graph(id='PCA-explained-variance', figure={}, config={'displayModeBar':False},
                                style={'fontSize': 14}),
                        #width=8, lg={'size': 8,  "offset": 0, 'order': 'first'}
                       # ),
                        ]),
        
        html.Br(),
        html.Br(),
        dbc.Row(dbc.Col(html.H1('Pathway analysis: Which genes and pathways cause the most variation in the data?',style={'textAlign':'center'}),
                        width={'size': 'true', 'offset': 0})),
        html.Hr(),
        dbc.Row(
                [
                    dbc.Col(dbc.Card(card_content_left, color="secondary", outline=True),),
                    #dbc.Col(dbc.Card(card_content_right, outline=False)),
                     ],),
        
        
        html.Hr(),
        
        dbc.Row(dbc.Col(dcc.Graph(id='Kegg-pathways-loadings', figure={}, config={'showTips':True}),
                       # width=8, lg={'size': 6,  "offset": 0, 'order': 'first'}
                        ),),
        dbc.Row(
            dbc.Col(
                dcc.Graph(id='cluster-heatmap', 
                #figure={'layout': go.Layout(margin={'t': 2})},
                config={'showTips':True}),
                        width={"offset": 0}),),
            ], className="bg-light")
# -------------------------------------------------------------------------------------
# Create PCA chart and bar chart for loadings
@app.callback(
    [Output(component_id='PCA-chart', component_property='figure'),
    Output(component_id='PCA-explained-variance', component_property='figure'),
    Output(component_id='Kegg-pathways-loadings', component_property='figure')],
    [Input(component_id='PCA_data_options', component_property='value'),
    Input(component_id='KEGG_missing_included',component_property='value'),
    Input(component_id='av_or_sum',component_property='value'),
    Input(component_id='loadings_for_graph',component_property='value')],
)

def PCA_for_graphs(PCA_data_options,KEGG_missing_included,av_or_sum,loadings_for_graph):
    PCA = countData.copy()
    PCA = make_df(PCA_data_options, PCA)
    if KEGG_missing_included == 'Remove "Missing" KEGG Pathways':
        PCA = PCA[PCA.KEGG_Pathway != 'Missing']

    df_genes = PCA.drop(labels={'Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL'}, axis=1)
    if len(df_genes.columns) < 5:
        to_plot_and_analyze, loadings, expl_var = PCA_script(df_genes, transpose = True, loadings = True, pcs = 4)
        to_plot_and_analyze = pd.DataFrame(to_plot_and_analyze, columns = ['PC1', 'PC2', 'PC3', 'PC4'])
        to_plot_and_analyze = pd.DataFrame(to_plot_and_analyze, columns = ['PC1', 'PC2', 'PC3', 'PC4'])
        pc = ['PC1','PC2','PC3','PC4']
        per_pc = ['PC','PC','PC','PC']
    else:
        to_plot_and_analyze, loadings, expl_var = PCA_script(df_genes, transpose = True, loadings = True, pcs = 5)
        to_plot_and_analyze = pd.DataFrame(to_plot_and_analyze, columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
        to_plot_and_analyze = pd.DataFrame(to_plot_and_analyze, columns = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])
        pc = ['PC1','PC2','PC3','PC4','PC5']
        per_pc = ['PC','PC','PC','PC','PC']

    PCA['loadings PC1'] = loadings[:,0]
    PCA['loadings PC2'] = loadings[:,1]
    PCA['loadings PC3'] = loadings[:,2]
    PCA['loadings PC4'] = loadings[:,3]
    to_plot_and_analyze['Condition'] = df_genes.columns.tolist()
    to_plot_and_analyze['Condition'] = to_plot_and_analyze['Condition'].str[:9]
    to_plot_and_analyze['Labels'] = df_genes.columns.tolist()
    to_plot_and_analyze['abs(PC1)'] = np.absolute(to_plot_and_analyze['PC2'])

    #chart PCA
    chart = px.scatter(to_plot_and_analyze, x="PC1", y="PC2", color="Condition",
                        hover_data={'Condition':False,'PC1':':.2f','PC2':':.2f','Labels':True},
                      #  width=800, height=500,
                        color_discrete_sequence=px.colors.qualitative.G10, template="simple_white",
                        )
    chart.update_yaxes(showgrid=True)
    chart.update_xaxes(showgrid=True)
    chart.update_layout(uniformtext_minsize=10, uniformtext_mode='hide')
    chart.update_traces(marker={'size': 15})

    # chart explained variance
    expl_var = pd.DataFrame(expl_var*100,columns=['Variance explained by each PC'])
    expl_var['PC'] = pc
    expl_var['percent PC'] = per_pc
    chart2 = px.bar(expl_var, x='percent PC',y='Variance explained by each PC',color='PC',
                    hover_data={'PC':False,'percent PC':False,'Variance explained by each PC':':.2f'},barmode='stack',width=300, height=450,
                    color_discrete_sequence=px.colors.qualitative.G10, template="simple_white",
                    )
    chart2.update_yaxes(visible=True,rangemode='tozero',showgrid=False)
    chart2.update_xaxes(visible=False,showgrid=True)


    #KEGG bar chart
    # choice_4 = st.radio("Loadings calculation", ('Average','Sum'))
    if av_or_sum == 'Average':
        PCA_KEGG_grouped = PCA.groupby('KEGG_Pathway').mean().reset_index()
    elif av_or_sum == 'Sum':
        PCA_KEGG_grouped = PCA.groupby('KEGG_Pathway').sum().reset_index()
    PCA_KEGG_grouped.melt()
    chart3 = px.bar(PCA_KEGG_grouped,y=loadings_for_graph,x='KEGG_Pathway',width=1250,height=650,
                    custom_data=['KEGG_Pathway'], color_discrete_sequence=px.colors.qualitative.G10,
                    template="simple_white",)
    chart3.update_layout(title='',xaxis={'categoryorder':'total descending'})
    chart3.update_xaxes(title_text='',side='top',tickangle=75,)
    chart3.update_yaxes(showgrid=True)
    #chart3.update_layout(hoverlabel=dict(bgcolor="white",font_size=12,font_family="Rockwell"))

    return chart, chart2, chart3


@app.callback(
    Output(component_id='cluster-heatmap', component_property='figure'),
    [Input(component_id='loadings_for_graph',component_property='value'),
    Input(component_id='KEGG_missing_included',component_property='value'),
    Input(component_id='PCA_data_options', component_property='value'),
    Input(component_id='Kegg-pathways-loadings', component_property='clickData')] #clickData is also good
)
def heatmap_PCA(loadings_for_graph,KEGG_missing_included,PCA_data_options,hover_dat):
    df_genes = countData.copy()
    df_genes = make_df(PCA_data_options, df_genes)
    if KEGG_missing_included == 'Remove "Missing" KEGG Pathways':
        df_genes = df_genes[df_genes.KEGG_Pathway != 'Missing']

    #print(f'hover data: {hover_dat}')
    #print(f'hover data: {filtered_kegg_pathway}')
    if hover_dat is None:
        y_labels = df_genes.Name
        df_genes = df_genes.drop(labels={'Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL'}, axis=1)
        to_plot_and_analyze, loadings, expl_var = PCA_script(df_genes, transpose = True, loadings = True, pcs = 5)
        df_genes['Name'] = y_labels.str[:28]
        df_genes['loadings PC1'] = loadings[:,0]
        df_genes['loadings PC2'] = loadings[:,1]
        df_genes['loadings PC3'] = loadings[:,2]
        df_genes['loadings PC4'] = loadings[:,3]

        heatmap = pd.DataFrame(df_genes.nlargest(40, loadings_for_graph))
        heatmap = heatmap.append(df_genes.nsmallest(40, loadings_for_graph))
        heatmap = heatmap.drop(labels={'loadings PC1','loadings PC2','loadings PC3','loadings PC4'}, axis=1)
        heatmap = heatmap.set_index('Name')
        columns = list(heatmap.columns.values)
        rows = list(heatmap.index)
        fig = clust(data=heatmap.transpose(), column_labels=rows,row_labels=columns,height=800,width=1300,
        color_threshold={
            'row': 250,
            'col': 700},
        color_list={
        'row': ['#636EFA', '#00CC96', '#19D3F3'],
        'col': ['#AB63FA', '#EF553B'],
        'bg': '#506784'},
        color_map= [
        [0.0, "#3366cc"],
        [0.5, "#FFFFFF"],
        [1.0, "#cb4827"]],
        display_ratio=[0.05, 0.15])
        fig.update_xaxes(tickangle=65,)
        fig.update_layout(margin=dict(l=10, r=20, t=0, b=20))
        return fig
    else:
        filtered_kegg_pathway = hover_dat['points'][0]['customdata']
        heatmap = df_genes[df_genes.KEGG_Pathway.isin(filtered_kegg_pathway)]
        heatmap = heatmap.drop(labels={'Fasta','ECNumber','KEGG_Pathway','KEGGSYMBOL'}, axis=1)
        heatmap = heatmap.set_index('Name')
        columns = list(heatmap.columns.values)
        rows = list(heatmap.index)
        fig = clust(data=heatmap.transpose(), column_labels=rows,row_labels=columns,height=800,width=1300,
        color_threshold={
            'row': 250,
            'col': 700},
        color_list={
        'row': ['#636EFA', '#00CC96', '#19D3F3'],
        'col': ['#AB63FA', '#EF553B'],
        'bg': '#506784'},
    color_map= [
        [0.0, "#3366cc"],
        [0.5, "#FFFFFF"],
        [1.0, "#cb4827"]],
        display_ratio=[0.05, 0.15])
        fig.update_xaxes(tickangle=65,)
        fig.update_layout(margin=dict(l=10, r=20, t=0, b=20))
        
        return fig


# make the data table, format it, and add to KEGG pathway page
    # could have where click on something and its gene information pops up.


# in the table, use 'clickData' and 'selectData' 
# to change output https://www.youtube.com/watch?v=G8r2BB3GFVY&list=RDCMUCqBFsuAz41sqWcFjZkqmJqQ&index=11 19 min



# somehow show the genome? -> just the fasta file maybe
    # the genome https://dash.gallery/dash-alignment-chart/
    # A gene table, and a reaction page



# add in all your comments as markdown at the bottom of each page, in an open-able tab
### 2) markdown_text = ''' ### Dash and Markdown Dash apps can be written in Markdown. Dash uses the [CommonMark](http://commonmark.org/) specification of Markdown. Check out their [60 Second Markdown Tutorial](http://commonmark.org/help/) if this is your first introduction to Markdown! ''' app.layout = html.Div([ dcc.Markdown(children=markdown_text) ]) ________

