from dash import dcc
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Input, Output, State
import plotly.express as px
import pandas as pd
import numpy as np
import pathlib
from app import app
from plotly.subplots import make_subplots
import plotly.graph_objects as go
from dash import dash_table
from cobra.flux_analysis import flux_variability_analysis

PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../datasets").resolve()

import cobra
import cobra.test
from cobra import Reaction, Metabolite
from cobra.io import read_sbml_model
from cobra.flux_analysis import production_envelope
from cobra.flux_analysis import flux_variability_analysis, pfba
from textwrap import shorten
from typing import List, Dict
from operator import attrgetter
from cobra.util.solver import linear_reaction_coefficients

#manipulating model
def set_yeast_extraction(model,ub=0,lb=0):
    # Block the uptake of 20 amino acids. Can also adjust uptake to whatever value want. Also block ammonia, lactate.
    amino_acids_inc = ['CYS_c','L-ASPARTATE_c','Glutamate','PHE_c','GLY_c','HIS_c','Isoleucine','LYS_c','Leucine','MET_c','ASN_c','Proline','Glutamine','ARG_c','SER_c','THR_c','VAL_c','AMMONIA','Tryptophan','L-ALPHA-ALANINE_c','TYR_c','L-LACTATE','Trehalose', 'ADENINE']
    for met_id in amino_acids_inc:
        exg_rxn = model.reactions.get_by_id('Exchange_'+met_id.replace('_c',''))
        exg_rxn.upper_bound = ub
        exg_rxn.lower_bound = lb
def test_PHA_production(model):
    pha_sink = cobra.Reaction('pha_sink')
    pha_sink.add_metabolites({model.metabolites.get_by_id('PHA_c'):-1})
    with model:
        model.add_reactions([pha_sink]) 
        model.objective = 'pha_sink'
        model.objective_direction = 'max'
        s1 = model.optimize()
    return s1
def test_ectoine_production(model):
    ect_sink = cobra.Reaction('ect_sink')
    ect_sink.add_metabolites({model.metabolites.get_by_id('ECTOINE_c'):-1})
    with model:
        model.add_reactions([ect_sink]) 
        model.objective = 'ect_sink'
        model.objective_direction = 'max'
        s1 = model.optimize()
    return s1
def testing_model(model):
    with model:
        solution=model.optimize()
        e=test_ectoine_production(model)
        p=test_PHA_production(model)
    return solution,e,p
def set_bound(model,rxn_id,lb=0,ub=0):
    rxn = model.reactions.get_by_id(rxn_id)
    rxn.lower_bound = lb
    rxn.upper_bound = ub
def reset_carbon(model,input_for_conditions,yeast_adjust):
    for c in carbon_sources:
        model.reactions.get_by_id(c).upper_bound = 0
    model.reactions.get_by_id(input_for_conditions).upper_bound = 10
    set_yeast_extraction(model,yeast_adjust,-yeast_adjust)

#functions for displaying fva for exchange reactions
def ddisplay_exchange_flux(_boundary_metabolites, frame, threshold=0.00):
    frame = frame.loc[
                (frame["Flux"].abs() >= threshold)
                | (frame["minimum"].abs() >= threshold)
                | (frame["maximum"].abs() >= threshold),
                :,
            ].copy()

    metabolites = {m.id: m for m in _boundary_metabolites}
    frame["Metabolite"] = [metabolites[met_id].name for met_id in frame["Metabolite"]]
    #frame["FVA Range"] = frame["minimum"].map('{:,.3f}'.format).astype(str)+" to "+frame["maximum"].map('{:,.3f}'.format).astype(str)
    return frame[["Metabolite","Flux",'minimum', 'maximum']]
def fva_exchange_summary(model, slider, threshold):
    coefficients = linear_reaction_coefficients(model)
    solution = pfba(model)
    fva = flux_variability_analysis(
                model=model,
                reaction_list=model.boundary,
                fraction_of_optimum=slider,
            )
    _objective: Dict["Reaction", float] = {
                rxn.copy(): coef for rxn, coef in coefficients.items()
            }
    _objective_value: float = sum(
                solution[rxn.id] * coef for rxn, coef in _objective.items()
            )
    _boundary: List["Reaction"] = [
            rxn.copy() for rxn in sorted(model.boundary, key=attrgetter("id"))
        ]
    _boundary_metabolites: List["Metabolite"] = [
            met.copy() for rxn in _boundary for met in rxn.metabolites
        ]
    flux = pd.DataFrame(
            data=[
                (rxn.id, met.id, rxn.get_coefficient(met.id), solution[rxn.id])
                for rxn, met in zip(_boundary, _boundary_metabolites)
            ],
            columns=["Rxn", "Metabolite", "factor", "Flux"],
            index=[r.id for r in _boundary],
        )
        # Scale fluxes by stoichiometric coefficient.
    flux["Flux"] *= flux["factor"]
    flux = flux.join(fva)
    view = flux[["Flux", "minimum", "maximum"]]
            # Set fluxes below model tolerance to zero.
    flux[["Flux", "minimum", "maximum"]] = view.where(
                view.abs() >= model.tolerance, 0
            )
            # Create the scaled compound flux.
    flux[["minimum", "maximum"]] = flux[["minimum", "maximum"]].mul(
                flux["factor"], axis=0
            )
            # Negative factors invert the minimum/maximum relationship.
    negative = flux["factor"] < 0
    tmp = flux.loc[negative, "maximum"]
    flux.loc[negative, "maximum"] = flux.loc[negative, "minimum"]
    flux.loc[negative, "minimum"] = tmp
            # Add zero to turn negative zero into positive zero for nicer display later.
    flux[["Flux", "minimum", "maximum"]] += 0
    # Create production table from producing fluxes or zero fluxes where the
        # metabolite is a product in the reaction.
    is_produced = (flux["Flux"] > 0)
    uptake_flux = flux.loc[
                is_produced, ["Flux", "minimum", "maximum", "Rxn", "Metabolite"]
            ].copy()
    is_consumed = (flux["Flux"] < 0)
    secretion_flux = flux.loc[
                is_consumed, ["Flux", "minimum", "maximum", "Rxn", "Metabolite"]
            ].copy()
    uptake_flux = ddisplay_exchange_flux(_boundary_metabolites, uptake_flux, threshold=0.00)
    secretion_flux = ddisplay_exchange_flux(_boundary_metabolites, secretion_flux, threshold=0.00)
    to_return = pd.concat([uptake_flux,secretion_flux])
    to_return['Optimal Flux'] = to_return['Flux'].map('{:,.4f}'.format)
    return to_return[["Metabolite","minimum", "maximum", "Optimal Flux"]].sort_values(by="Optimal Flux", ascending=False)

#functions for displaying fva for metabolites
def ddisplay_flux(_reactions, frame, threshold= 0.00):

        frame = frame.loc[
                (frame["Flux"].abs() >= threshold)
                | (frame["minimum"].abs() >= threshold)
                | (frame["maximum"].abs() >= threshold),
                :,
            ].copy()
        reactions = {r.id: r for r in _reactions}
        frame["Stoichiometry"] = [
            reactions[rxn_id].build_reaction_string(False)
            for rxn_id in frame["Rxn"]
        ]
        frame["FVA Range"] = "("+frame["minimum"].map('{:,.3f}'.format).astype(str)+" - "+frame["maximum"].map('{:,.3f}'.format).astype(str)+"), optimal: "+frame["Flux"].map('{:,.4f}'.format).astype(str)

        to_return = frame[["%", "FVA Range","Rxn", "Stoichiometry"]]
        to_return['%'] = to_return['%'].map('{:,.2f}'.format)
        #to_return['Optimal Flux'] = to_return['Flux'].map('{:,.4f}'.format)
        return to_return[["%","FVA Range","Rxn", "Stoichiometry"]]
def metabolite_flux_summary(model, metabolite_summary, thresh, slider):
    # https://cobrapy.readthedocs.io/en/latest/_modules/cobra/summary/metabolite_summary.html#MetaboliteSummary
        _reactions: List["Reaction"] = [
                r.copy() for r in sorted(model.metabolites.get_by_id(metabolite_summary).reactions, key=attrgetter("id"))
            ]
        solution = pfba(model)
        fva = flux_variability_analysis(
                model=model,
                reaction_list=[r.id for r in _reactions],
                fraction_of_optimum=slider,
            )
        flux = pd.DataFrame(
                data=[
                    (r.id,solution[r.id],r.get_coefficient(metabolite_summary),)
                    for r in _reactions],
                columns=["Rxn", "Flux", "factor"],
                index=[r.id for r in _reactions],
            )
        flux["Flux"] *= flux["factor"]
        
        flux = flux.join(fva)
        view = flux[["Flux", "minimum", "maximum"]]
            # Set fluxes below model tolerance to zero.
        flux[["Flux", "minimum", "maximum"]] = view.where(
                view.abs() >= model.tolerance, 0
            )
            # Create the scaled compound flux.
        flux[["minimum", "maximum"]] = flux[["minimum", "maximum"]].mul(
                flux["factor"], axis=0
            )
            # Negative factors invert the minimum/maximum relationship.
        negative = flux["factor"] < 0
        tmp = flux.loc[negative, "maximum"]
        flux.loc[negative, "maximum"] = flux.loc[negative, "minimum"]
        flux.loc[negative, "minimum"] = tmp
            # Add zero to turn negative zero into positive zero for nicer display later.
        flux[["Flux", "minimum", "maximum"]] += 0

        # Create production table from producing fluxes or zero fluxes where the metabolite is a product in the reaction.
        if thresh == 0:
            is_produced = (flux["Flux"] > 0) | ((flux["Flux"] == 0) & (flux["factor"] > 0))
        else:
            is_produced = (flux["Flux"] > 0)
        producing_flux = flux.loc[
                is_produced, ["Flux", "minimum", "maximum", "Rxn"]
            ].copy()
        production = producing_flux["Flux"].abs()
        producing_flux["%"] = production / production.sum()
        # Create consumption table from consuming fluxes or zero fluxes where the metabolite is a substrate in the reaction.
        if thresh == 0:
            is_consumed = (flux["Flux"] < 0) | ((flux["Flux"] == 0) & (flux["factor"] < 0))
        else:
            is_consumed = (flux["Flux"] < 0)        
        consuming_flux = flux.loc[
                is_consumed, ["Flux", "minimum", "maximum", "Rxn"]
            ].copy()
        consumption = consuming_flux["Flux"].abs()
        consuming_flux["%"] = consumption / consumption.sum()

        producing_flux = ddisplay_flux(_reactions, producing_flux, 0.00)
        consuming_flux = ddisplay_flux(_reactions, consuming_flux, 0.00)
        producing_flux = producing_flux.rename(columns={"Stoichiometry": "Producing Rxn"})
        consuming_flux = consuming_flux.rename(columns={"Stoichiometry": "Consuming Rxn"})

        return producing_flux.sort_values(by="%",ascending=False), consuming_flux.sort_values(by="%",ascending=False)

#Charts 
def make_chart_varying_media(model,glcs,media):
    with model:
        rgs_s = []
        rgs_e = []
        rgs_p = []
        for glc in glcs:
            model.reactions.get_by_id(media).lower_bound = glc
            model.reactions.get_by_id(media).upper_bound = glc
            s3 = model.optimize()
            e3 = test_ectoine_production(model)
            p3 = test_PHA_production(model)
            if s3.objective_value is None:
                rgs_s.append(0)
            else:
                rgs_s.append(s3.objective_value)
            if e3.objective_value is None:
                rgs_e.append(0)
            else:
                rgs_e.append(e3.objective_value)
            if p3.objective_value is None:
                rgs_p.append(0)
            else:
                rgs_p.append(p3.objective_value)
        return rgs_s, rgs_e, rgs_p
def figure_for_production(model,carbon_source,objective_name):
    lb=0
    ub=20
    a=1
    df = pd.DataFrame(columns={carbon_source,'Growth','Ectoine','PHA'})
    df[carbon_source] = np.arange(lb,ub,a)
    df.Growth, df.Ectoine, df.PHA = make_chart_varying_media(model,df[carbon_source],carbon_source)

    # Create figure with secondary y-axis
    fig = make_subplots(subplot_titles=("Carbon Uptake and Cell Production and Growth", "TD1.0 Production and Growth"), specs=[[{"secondary_y": True}]])
    # Add traces
    fig.add_trace(go.Scatter(x=df[carbon_source],y=df.Growth,name="Objective"), secondary_y=False)
    fig.add_trace(go.Scatter(x=df[carbon_source],y=df.Ectoine,name="Ectoine"), secondary_y=True,)
    fig.add_trace(go.Scatter(x=df[carbon_source],y=df.PHA,name="PHA"), secondary_y=True)

    fig.update_layout(template="simple_white", )
    fig.update_xaxes(title_text=carbon_source)
    fig.update_yaxes(title_text=objective_name, secondary_y=False)
    fig.update_yaxes(title_text="Ectoine, PHA Production/CDW", secondary_y=True)
    return fig
def figure_for_flux_exchange(df):
    fig = go.Figure(data=[go.Candlestick(
                x=df['Metabolite'],open=df['Optimal Flux'], high=df['maximum'],low=df['minimum'], close=df['Optimal Flux'])])    
    fig.update_layout(xaxis_rangeslider_visible=False)
    fig.update_yaxes(title_text="Flux")
    fig.update_layout(template="simple_white", title="Consumption and Production Reactions (FVA)", title_x=0.5)
    return fig

def phenotypic_phase_plane(model, var1, var2, carbon_source ,objective, number_1, number_2, number_3, number_4):
    with model:
        if (number_1 is not None and number_2 is not None and number_3 is not None and number_4 is not None):
            set_bound(model, var1,lb=number_1, ub=number_3)
            set_bound(model, var2,lb=number_2, ub=number_4)
        
        prod_env = production_envelope(
            model, [var2, var1], objective = objective, carbon_sources=carbon_source)

    Zmask = np.isfinite(prod_env['flux_maximum'])
    prod_env = prod_env[[var2, var1,'flux_maximum']]
    prod_env = prod_env[Zmask]
    fig = go.Figure(data=[go.Mesh3d(z=prod_env['flux_maximum'], x=prod_env[var2], y=prod_env[var1],
                                opacity=0.8,colorscale=[[0, 'gold'],[0.5, 'orange'],[1, 'red']],
                                    colorbar_title='flux',intensity=prod_env['flux_maximum'],showscale=False)])

    fig.update_layout(title='Phenotypic Phase Plane',title_x=0.5, scene = {
                   "xaxis_title_text":var2,"yaxis_title_text":var1,"zaxis_title_text":objective,
                "aspectratio": {"x": 1, "y": 1, "z": 0.6}}, 
                  template="plotly_white",width=700, height=700,
                 #margin=dict(l=10, r=10, b=100, t=2),
                 )
                 #scene_camera_eye=dict(x=0.88, y=-0.88, z=0.14))
    fig.update_traces(contour=dict(show=True, color="limegreen", ))
    return fig


carbon_sources = ['Exchange_Glucopyranose','Exchange_L-LACTATE','Exchange_ACET','Exchange_CIT','GLYCtex','Exchange_Trehalose','BUTtex','EX_SUC_e','ETOHtex', 'FRUtex']
model2 = read_sbml_model(str(DATA_PATH.joinpath('Gap_filled_media_chalmers_v2.sbml')))
info_content_right = html.Div([
    dbc.Button( "  Guide  ", style={'fontSize': 14, "horizontalAlign": "middle", 'margin':'auto', 'textAlign':'center'},id="hover-target-2"),
            dbc.Modal(
            [
                dbc.ModalHeader(dbc.ModalTitle("Halomonas TD1.0 GEM model Guide",), close_button=True,),
                dbc.ModalBody(dcc.Markdown('''
                ### It can be a little slow to load; If it takes longer than 30 seconds, switch to LB+ media, then switch back!
                ###  
                #### Change the carbon source, media, objective function, FVA ratio, and introduce any KOs.
                - To **refresh the model**, refresh the page.
                - The **objective function** is the reaction the model is optimizing for. By default this is Biomass, but you can optimize for any product.
                    - For the TD1.0 Production and Growth graph, the objective reaction of choice is the **blue line**. The Ectoine and PHA lines are solutions where ectoine and pha are the objective value, respectively.
                - At this time you can **introduce KO** on a gene basis into the model. A gene KO will also KO any reaction in the model that it is necessary for.
                - LB+ Media allows the consumption/production of all 20 amino acids (Not recommended unless interested in AA consumption.)
                ####  
                #### The model rarely has one single solution. FVA (flux variability analysis) allows us to explore the solution space. FVA will:
                - Solve for the optimum objective (using parsimonious FBA, pFBA)
                - Set the objective value to be at the ratio specified. For example, we set biomass to be 95% of optimum value.
                - Next, the model **iterates to find the minimum and maximum flux** in the range of possible solutions.
                - Theoretically, reactions with small variation tend to be more important.
                ####  
            ''', style={'fontSize': 14}), ),
            
            ],
            id="modal-centered", size="lg",
            centered=True,
            is_open=False,   
        ),])
button_two = html.Div([
    dbc.Button(
            "Info", style={'fontSize': 14},
            id="hover-target-1", 
            className="me-1",
            n_clicks=0,
        ),
        dbc.Popover(
            dcc.Markdown('''
            - Hit **enter** after final entry. Graph may not populate if default reaction bounds are (0, 0)
            - May take ~20 seconds to load as it is computationally taxing!
            ''', 
            style={'fontSize': 12}), target="hover-target-1",placement="bottom",offset="50,20",body=True,trigger="hover", 
        ),

])

layout = html.Div([
        html.Br(),
        dbc.Row(dbc.Col(html.H1("Halomonas GEM Model: Flux Analysis", style={'textAlign':'center','fontSize': 30}))),
        html.Br(),
        html.Hr(),
        dbc.Row([dbc.Col(html.H4("Select carbon source and media:", style={'fontSize': 16}),
                        width={'size': 3, 'offset': 1}),
                dbc.Col(html.H4("Select reaction for objective function:", style={'fontSize': 16}),
                        width={'size': 4, 'offset': 0}),
                dbc.Col(html.H4("Select gene to KO of model:", style={'fontSize': 16}),
                        width={'size': 4, 'offset': 0})
                        ], className="mx-2"),
          
        dbc.Row([
            dbc.Col(info_content_right,width={'size': 1, "offset": 0},),
            dbc.Col(dcc.Dropdown(id='carbon_source', value='Exchange_Glucopyranose',clearable=False,
                                    persistence=True, persistence_type='memory',
                                    options=[{'label': x, 'value': x} for x in carbon_sources],
                                    style={'fontSize': 14}), width={'size': 2, "offset": 0}),
            dbc.Col(dbc.RadioItems(id='Media_selection',value=0,
                                    options=[
                                        {'label':'Minimal_Media   ','value':0},
                                        {'label':'LB+ Media','value':10}],
                                    style={'fontSize': 12}),
                                    width={'size': 1, "offset": 0},),
            
            dbc.Col(dcc.Dropdown(id='objective', value='Biomass_v1',clearable=False,
                                    options=[{'label': x.name, 'value': x.id} for x in model2.reactions],
                                    style={'fontSize': 14}), width={'size': 4}),
            dbc.Col(dcc.Dropdown(id='knock_out',value='Nothing',clearable=False,
                                    options=[{'label': x.name, 'value': x.id} for x in model2.genes],
                                    style={'fontSize': 14}), width={'size': 4}),
            
 
            ], className="mx-2"),
            html.Hr(),
            
            #graphs
            dbc.Row([
                dbc.Col([
                    dcc.Graph(id='growth', figure={},config={'showTips':True}), html.Hr(),
                    dcc.Graph(id='model_summary',  figure={},config={'showTips':True}), html.Hr(),        
                    dbc.Row([dbc.Col(html.H4("Select variables for flux envelope:", style={'fontSize': 16}),
                                width={'size': 6, 'offset': 0}),
                            dbc.Col(html.H4("Adjust Bounds for reactions (hit enter):", style={'fontSize': 16}),
                                width={'size': 5, 'offset': 0}),
                                ], className="mx-2"),
                    dbc.Row([dbc.Col([dcc.Dropdown(id='variable_1', value='Exchange_OXYGEN-MOLECULE',clearable=False,
                                        options=[{'label': x.name, 'value': x.id} for x in model2.reactions],style={'fontSize': 14}),
                                    dcc.Dropdown(id='variable_2', value='Exchange_Glucopyranose', clearable=False,
                                        options=[{'label': x.name, 'value': x.id} for x in model2.reactions],style={'fontSize': 14}),
                                        ], width={'size': 6}),
                            dbc.Col([dcc.Input(id="number_1",type='number',debounce=True, placeholder="min",size='10',min=0, max=100, step=1),
                                    dcc.Input(id="number_2",type='number',debounce=True, placeholder="min",size='10',min=0, max=100, step=1),
                            ],  width={'size': 1,'offset':0}),
                            dbc.Col([dcc.Input(id="number_3",type='number',debounce=True, placeholder="max",size='10',min=0, max=100, step=1),
                                    dcc.Input(id="number_4",type='number',debounce=True, placeholder="max",size='10',min=0, max=100, step=1),
                            ],  width={'size': 1}),
                            dbc.Col(button_two,width={'size': 1, "offset": 1}),
                    ], className="mx-2"),
                    html.Hr(),
                    dcc.Graph(id='pheno_plane', figure={}, config={'showTips':True}),], width={'size': 7, "offset": 0}),
                    
                    
            dbc.Col([html.Br(), html.Hr(),
                        dbc.Row(dbc.Col(html.H4("Select Metabolite for consuming and producing reactions.", style={'fontSize': 16}),), className="mx-2"),
                        dbc.Row(dbc.Col(dcc.Dropdown(id='metabolite_summary', value='ATP_c',clearable=False,
                                    persistence=True, persistence_type='memory',
                                    options=[{'label': x.name, 'value': x.id} for x in model2.metabolites],
                                    style={'fontSize': 14}))),html.Br(),
                        dbc.Row([dbc.Col(dbc.RadioItems(id='threshold',value=0.0000000000001,
                                    options=[{'label':'All Reactions','value':0},{'label':'Only Non-Zero Reactions','value':0.0000000000001}],
                                    style={'fontSize': 14}), width={'size': 5}),   
                                dbc.Col([html.H4("Flux Variability % optimal solution", style={'text-align':'center','fontSize': 14}),
                                    dcc.Slider(min=0.1, max=1, step=0.05, value=0.95, marks = None, id='my-slider',
                                        tooltip={"placement": "bottom", "always_visible": True})], width={'size': 7}) ,      
                                ]),
                    html.Hr(),
                    dbc.Col(dash_table.DataTable(id='met_prod_summary', 
                        page_size=15, filter_action="native", sort_action="native",
                        style_header={'backgroundColor': 'rgb(247, 211, 207)','color': 'black','fontWeight': 'bold','fontSize': 14, 'fontType':'Times'},
                        style_cell_conditional=[ {'textAlign': 'left', 'fontSize': 12, 'fontType':'Times'}],
                        style_data={  'whiteSpace': 'normal', 'height': 'auto'}, css=[{"selector": "table", "rule": "width: 100%;"}]), 
                            #width={'size': 4, "offset": 0}
                            ),
                    html.Hr(),
                    dbc.Col(dash_table.DataTable(id='met_con_summary', 
                        page_size=15, filter_action="native", sort_action="native",
                        style_header={'backgroundColor': 'rgb(247, 211, 207)','color': 'black','fontWeight': 'bold','fontSize': 14, 'fontType':'Times'},
                        style_cell_conditional=[ {'textAlign': 'left', 'fontSize': 12, 'fontType':'Times'}], 
                        style_data={ 'whiteSpace': 'normal', 'height': 'auto'
                                        }, css=[{"selector": "table", "rule": "width: 100%;"}]), 
                            #width={'size': 4, "offset": 0}
                    )], className="m-8", width={'size': 5, "offset": 0}),
            ]),
        ])

@app.callback(Output("modal-centered", "is_open"),
            Input("hover-target-2", "n_clicks"),
            State("modal-centered", "is_open"))
def info_card(n1, is_open):
    if n1:
        return not is_open
    return is_open


@app.callback(
     [Output('growth', 'figure'),
     Output('model_summary', 'figure'),
     Output('pheno_plane', 'figure'),
     Output(component_id='met_prod_summary', component_property='data'), Output(component_id='met_prod_summary', component_property='columns'),
     Output(component_id='met_con_summary', component_property='data'), Output(component_id='met_con_summary', component_property='columns')],
     [Input('carbon_source', 'value'),
     Input('Media_selection','value'),
     Input('objective','value'),
     Input('knock_out','value'),
     Input('metabolite_summary', 'value'),
     Input('threshold', 'value'),
     Input('my-slider', 'value'),
     Input('variable_1','value'), Input('variable_2','value'), Input('number_1','value'), Input('number_2','value'), Input('number_3','value'), Input('number_4','value')])
def flux_the_gem(carbon_source, Media_selection, objective, knock_out, metabolite_summary, threshold, slider, 
                        plane_var1, plane_var2, number_1, number_2, number_3, number_4):
    model = model2.copy()
    
    #reset the carbon source
    set_bound(model,'EX_HEXANOATE_e',0,0)
    set_bound(model,'Exchange_UREA',lb=-1,ub=10)
    reset_carbon(model,carbon_source,Media_selection)
    
    #change objective when user changes
    model.objective = objective
    # The upper bound should be 1000, so that we get  the actual optimal value
    model.reactions.get_by_id(objective).upper_bound = 1000.
    
    #KOs for model
    if knock_out != "Nothing" and knock_out != None:
        model.genes.get_by_id(knock_out).knock_out()
    #growth, ectoine, pha = testing_model(model)    

    #exchange reactions
    df = fva_exchange_summary(model, slider, threshold)

    # metabolite summary
    producing, consuming = metabolite_flux_summary(model, metabolite_summary, threshold, slider)
    produce_column, produce_data = [{'name': col, 'id': col} for col in producing.columns], producing.to_dict(orient='records')
    consuming_column, consuming_data = [{'name': col, 'id': col} for col in consuming.columns], consuming.to_dict(orient='records')

    return (figure_for_production(model,carbon_source, objective), figure_for_flux_exchange(df), 
            phenotypic_phase_plane(model, plane_var1,plane_var2, carbon_source, objective, number_1, number_2, number_3, number_4),
                    produce_data, produce_column, 
                    consuming_data, consuming_column)



    