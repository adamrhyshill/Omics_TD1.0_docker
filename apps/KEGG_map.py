from dash import dcc
from dash import html
import dash_bootstrap_components as dbc
from dash.dependencies import Input, Output, State
import plotly.express as px
#import plotly.io as pio
import numpy as np
import pandas as pd
import pathlib
from app import app
from Bio.KEGG.REST import kegg_get, kegg_list
from Bio.KEGG.KGML import KGML_parser
import os
import io
import base64

#-------------
import tempfile
from reportlab.lib import colors
from reportlab.pdfgen import canvas
from PIL import Image
from io import BytesIO
from urllib.request import urlopen

#------------needed functions for KGML
def darken(color, factor=0.7):
    """Return darkened color as a ReportLab RGB color.
    Take a passed color and returns a Reportlab color that is darker by the
    factor indicated in the parameter.
    """
    newcol = color_to_reportlab(color)
    for a in ["red", "green", "blue"]:
        setattr(newcol, a, factor * getattr(newcol, a))
    return newcol
def color_to_reportlab(color):
    """Return the passed color in Reportlab Color format.
    We allow colors to be specified as hex values, tuples, or Reportlab Color
    objects, and with or without an alpha channel. This function acts as a
    Rosetta stone for conversion of those formats to a Reportlab Color
    object, with alpha value.
    Any other color specification is returned directly
    """
    # Reportlab Color objects are in the format we want already
    if isinstance(color, colors.Color):
        return color
    elif isinstance(color, str):  # String implies hex color
        if color.startswith("0x"):  # Standardise to octothorpe
            color.replace("0x", "#")
        if len(color) == 7:
            return colors.HexColor(color)
        else:
            try:
                return colors.HexColor(color, hasAlpha=True)
            except TypeError:  # Catch pre-2.7 Reportlab
                raise RuntimeError(
                    "Your reportlab seems to be too old, try 2.7 onwards"
                ) from None
    elif isinstance(color, tuple):  # Tuple implies RGB(alpha) tuple
        return colors.Color(*color)
    return color
def get_temp_imagefilename(url):
    """Return filename of temporary file containing downloaded image.
    Create a new temporary file to hold the image file at the passed URL
    and return the filename.
    """
    img = urlopen(url).read()
    im = Image.open(BytesIO(img))
    f = tempfile.NamedTemporaryFile(delete=False, suffix=".png")
    fname = f.name
    f.close()
    im.save(fname, "PNG")
    return fname

#-------------------- Class for KGML
class KGMLCanvas:
    """Reportlab Canvas-based representation of a KGML pathway map."""

    def __init__(
        self,
        pathway,
        import_imagemap=False,
        label_compounds=True,
        label_orthologs=True,
        label_reaction_entries=True,
        label_maps=True,
        show_maps=False,
        fontname="Helvetica",
        fontsize=6,
        draw_relations=True,
        show_orthologs=True,
        show_compounds=True,
        show_genes=True,
        show_reaction_entries=True,
        margins=(0.02, 0.02),
    ):
        """Initialize the class."""
        self.pathway = pathway
        self.show_maps = show_maps
        self.show_orthologs = show_orthologs
        self.show_compounds = show_compounds
        self.show_genes = show_genes
        self.show_reaction_entries = show_reaction_entries
        self.label_compounds = label_compounds
        self.label_orthologs = label_orthologs
        self.label_reaction_entries = label_reaction_entries
        self.label_maps = label_maps
        self.fontname = fontname
        self.fontsize = fontsize
        self.draw_relations = draw_relations
        self.non_reactant_transparency = 0.3
        self.import_imagemap = import_imagemap  
        self.margins = margins

    def draw(self, filename):
        """Add the map elements to the drawing."""
        imfilename = get_temp_imagefilename(self.pathway.image)
        im = Image.open(imfilename)
        cwidth, cheight = im.size
        # Instantiate canvas
        self.drawing = canvas.Canvas(filename,bottomup=0,pagesize=(cwidth * (1 + 2 * self.margins[0]),cheight * (1 + 2 * self.margins[1]),),)
        self.drawing.setFont(self.fontname, self.fontsize)
        # Transform the canvas to add the margins
        self.drawing.translate(self.margins[0] * self.pathway.bounds[1][0],self.margins[1] * self.pathway.bounds[1][1],)
        # Add the map image, if required
        if self.import_imagemap:
            self.drawing.saveState()
            self.drawing.scale(1, -1)
            self.drawing.translate(0, -cheight)
            self.drawing.drawImage(imfilename, 0, 0)
            self.drawing.restoreState()
        # Add the reactions, compounds and maps
        # Maps go on first, to be overlaid by more information.
        # By default, they're slightly transparent.
        if self.show_maps:
            self.__add_maps()
        if self.show_reaction_entries:
            self.__add_reaction_entries()
        if self.show_orthologs:
            self.__add_orthologs()
        if self.show_compounds:
            self.__add_compounds()
        if self.show_genes:
            self.__add_genes()
        # TODO: complete draw_relations code
        # if self.draw_relations:
        #    self.__add_relations()
        # Write the pathway map to PDF
        self.drawing.save()

    def __add_maps(self):
        """Add maps to the drawing of the map (PRIVATE).
        We do this first, as they're regional labels to be overlaid by
        information.  Also, we want to set the color to something subtle.
        We're using Hex colors because that's what KGML uses, and
        Reportlab doesn't mind.
        """
        for m in self.pathway.maps:
            for g in m.graphics:
                self.drawing.setStrokeColor("#888888")
                self.drawing.setFillColor("#DDDDDD")
                self.__add_graphics(g)
                if self.label_maps:
                    self.drawing.setFillColor("#888888")
                    self.__add_labels(g)
    def __add_graphics(self, graphics):
        """Add the passed graphics object to the map (PRIVATE).
        Add text, add after the graphics object, for sane Z-ordering.
        """
        if graphics.type == "line":
            p = self.drawing.beginPath()
            x, y = graphics.coords[0]
            # There are optional settings for lines that aren't necessarily
            # part of the KGML DTD
            if graphics.width is not None:
                self.drawing.setLineWidth(graphics.width)
            else:
                self.drawing.setLineWidth(1)
            p.moveTo(x, y)
            for (x, y) in graphics.coords:
                p.lineTo(x, y)
            self.drawing.drawPath(p)
            self.drawing.setLineWidth(1)  # Return to default
        # KGML defines the (x, y) coordinates as the centre of the circle/
        # rectangle/roundrectangle, but Reportlab uses the coordinates of the
        # lower-left corner for rectangle/elif.
        if graphics.type == "circle":
            self.drawing.circle(
                graphics.x, graphics.y, graphics.width * 0.5, stroke=1, fill=1
            )
        elif graphics.type == "roundrectangle":
            self.drawing.roundRect(
                graphics.x - graphics.width * 0.5,
                graphics.y - graphics.height * 0.5,
                graphics.width,
                graphics.height,
                min(graphics.width, graphics.height) * 0.1,
                stroke=1,
                fill=1,
            )
        elif graphics.type == "rectangle":
            self.drawing.rect(
                graphics.x - graphics.width * 0.5,
                graphics.y - graphics.height * 0.5,
                graphics.width,
                graphics.height,
                stroke=1,
                fill=1,
            )
    def __add_labels(self, graphics):
        """Add labels for the passed graphics objects to the map (PRIVATE).
        We don't check that the labels fit inside objects such as circles/
        rectangles/roundrectangles.
        """
        if graphics.type == "line":
            # We use the midpoint of the line - sort of - we take the median
            # line segment (list-wise, not in terms of length), and use the
            # midpoint of that line.  We could have other options here,
            # maybe even parameterising it to a proportion of the total line
            # length.
            mid_idx = len(graphics.coords) * 0.5
            if not int(mid_idx) == mid_idx:
                idx1, idx2 = int(mid_idx - 0.5), int(mid_idx + 0.5)
            else:
                idx1, idx2 = int(mid_idx - 1), int(mid_idx)
            x1, y1 = graphics.coords[idx1]
            x2, y2 = graphics.coords[idx2]
            x, y = 0.5 * (x1 + x2), 0.5 * (y1 + y2)
        elif graphics.type == "circle":
            x, y = graphics.x, graphics.y
        elif graphics.type in ("rectangle", "roundrectangle"):
            x, y = graphics.x, graphics.y
        # How big so we want the text, and how many characters?
        if graphics._parent.type == "map":
            text = graphics.name
            self.drawing.setFont(self.fontname, self.fontsize + 2)
        elif len(graphics.name) < 15:
            text = graphics.name
        else:
            text = graphics.name[:12] + "..."
        self.drawing.drawCentredString(x, y, text)
        self.drawing.setFont(self.fontname, self.fontsize)
    def __add_orthologs(self):
        """Add 'ortholog' Entry elements to the drawing of the map (PRIVATE).
        In KGML, these are typically line objects, so we render them
        before the compound circles to cover the unsightly ends/junctions.
        """
        for ortholog in self.pathway.orthologs:
            for g in ortholog.graphics:
                self.drawing.setStrokeColor(color_to_reportlab(g.fgcolor))
                self.drawing.setFillColor(color_to_reportlab(g.bgcolor))
                self.__add_graphics(g)
                if self.label_orthologs:
                    # We want the label color to be slightly darker
                    # (where possible), so it can be read
                    self.drawing.setFillColor(darken(g.fgcolor))
                    self.__add_labels(g)
    def __add_reaction_entries(self):
        """Add Entry elements for Reactions to the map drawing (PRIVATE).
        In KGML, these are typically line objects, so we render them
        before the compound circles to cover the unsightly ends/junctions
        """
        for reaction in self.pathway.reaction_entries:
            for g in reaction.graphics:
                self.drawing.setStrokeColor(color_to_reportlab(g.fgcolor))
                self.drawing.setFillColor(color_to_reportlab(g.bgcolor))
                self.__add_graphics(g)
                if self.label_reaction_entries:
                    # We want the label color to be slightly darker
                    # (where possible), so it can be read
                    self.drawing.setFillColor(darken(g.fgcolor))
                    self.__add_labels(g)
    def __add_compounds(self):
        """Add compound elements to the drawing of the map (PRIVATE)."""
        for compound in self.pathway.compounds:
            for g in compound.graphics:
                # Modify transparency of compounds that don't participate
                # in reactions
                fillcolor = color_to_reportlab(g.bgcolor)
                if not compound.is_reactant:
                    fillcolor.alpha *= self.non_reactant_transparency
                self.drawing.setStrokeColor(color_to_reportlab(g.fgcolor))
                self.drawing.setFillColor(fillcolor)
                self.__add_graphics(g)
                if self.label_compounds:
                    if not compound.is_reactant:
                        t = 0.3
                    else:
                        t = 1
                    self.drawing.setFillColor(colors.Color(0.2, 0.2, 0.2, t))
                    self.__add_labels(g)
    def __add_genes(self):
        """Add gene elements to the drawing of the map (PRIVATE)."""
        for gene in self.pathway.genes:
            for g in gene.graphics:
                self.drawing.setStrokeColor(color_to_reportlab(g.fgcolor))
                self.drawing.setFillColor(color_to_reportlab(g.bgcolor))
                self.__add_graphics(g)
                if self.label_compounds:
                    self.drawing.setFillColor(darken(g.fgcolor))
                    self.__add_labels(g)
#----------------------------------- Functions for PCA

def standardize_2(df):
    df = df[df>0]
    df = np.log2(df)
    df = (df - df.mean())/df.std()
    return df

def make_df(choice, PCA):
    if choice == 'all_transcript':
        PCA = PCA[['Fasta','Name','ECNumber','KEGGSYMBOL','60gLNaCL_RPKM', '20gLNaCL_RPKM','100gLNaCL_RPKM','Urea3.6gL_RPKM','09Fermentation_RPKM','19Fermentation_RPKM', '30Fermentation_RPKM']].copy()
    elif choice == 'all_protein':
        PCA = PCA[['Fasta','Name','ECNumber','KEGGSYMBOL','60gLNaCL_ProteinpgDW','20gLNaCL_ProteinpgDW','100gLNaCL_ProteinpgDW','Urea3.6gL_ProteinpgDW', '09Fermentation_ProteinpgDW', '19Fermentation_ProteinpgDW','30Fermentation_ProteinpgDW']].copy()        
    for column in PCA.columns.drop(['Fasta','Name','ECNumber','KEGGSYMBOL']):
        PCA[column] = standardize_2(PCA[column])
    PCA = PCA.dropna()
    return PCA

def its_a_kegg_pathway(path_input):
    pathstring = kegg_get(path_input).read()
    try:
        pathstring = pathstring[:pathstring.index("REFERENCE")]
    except ValueError:
        pathstring = pathstring
    pathstring = pathstring[pathstring.index('ORTHOLOGY') + len('ORTHOLOGY'):]
    
    try:
        pathstring_compound = pathstring[:pathstring.index("COMPOUND")]
        pathstring_compound = ' - ' + pathstring_compound.replace("\n            ", "\n - ")
        pathstring_gene = pathstring[pathstring.index('COMPOUND') + len('COMPOUND'):]
        pathstring_gene = ' - ' + pathstring_gene.replace("\n            ", "\n - ")
    except:
        pathstring_compound = "No Genes"
    return pathstring_compound, (' - ' + pathstring.replace("\n            ", "\n - "))

PATH = pathlib.Path(__file__).parent
DATA_PATH = PATH.joinpath("../datasets").resolve()

countData = pd.read_csv(DATA_PATH.joinpath("R_Landing.csv"))  
t_or_p = ['Transcriptomics','Proteomics']
PCA_data_options = ['LB60','LB20','LB100','HighN','9hr','19hr','30hr','30 - 9hr', 'LB20 - LB60', 'LB100 - LB60']

# ------------------------------------------------------------
pathways = pd.read_table(io.StringIO((kegg_list('pathway', 'ko').read())),names=['KEGG_Gene','Name'])
outdir = '/Users/hellpark/Desktop/'
# kegg 
# https://github.com/biopython/biopython/blob/master/Bio/KEGG/REST.py 

card_content_compound = [dbc.CardHeader("Genes", style={'fontSize': 16}),
    dbc.CardBody([html.P(html.Div(id='compound-for-kegg', style={'whiteSpace': 'pre-line', 'fontSize': 14,"maxHeight": "350px", "overflow": "scroll"}),
                className="card-text",),]),]
card_content_gene = [dbc.CardHeader("Compounds",style={'fontSize': 16}),
    dbc.CardBody([html.P(html.Div(id='gene-for-kegg', style={'whiteSpace': 'pre-line', 'fontSize': 14,"maxHeight": "350px", "overflow": "scroll"}),
                className="card-text",)],),]
card_content_colors = [
            dbc.Col(html.Div("Info", style={'color': '#060EFB','fontSize': 20,'textAlign':'center', 'background-color': '#060EFB','t': 0},className="m-0")),
            dbc.Col(html.Div("Info", style={'color': '#4449F0','fontSize': 20,'textAlign':'center', 'background-color': '#4449F0'},className="m-0")),
            dbc.Col(html.Div("Info", style={'color': '#B2AAF1','fontSize': 20,'textAlign':'center', 'background-color': '#B2AAF1'})),
            dbc.Col(html.Div("Info", style={'color': '#c4fdf0','fontSize': 20,'textAlign':'center', 'background-color': '#c4fdf0'})),
            dbc.Col(html.Div("Info", style={'color': '#fbfb8d','fontSize': 20,'textAlign':'center', 'background-color': '#fbfb8d'})),
            dbc.Col(html.Div("Info", style={'color': '#FDA88C','fontSize': 20,'textAlign':'center', 'background-color': '#FDA88C'})),
            dbc.Col(html.Div("Info", style={'color': '#FC5749','fontSize': 20,'textAlign':'center', 'background-color': '#FC5749'})),
            dbc.Col(html.Div("Info", style={'color': '#FA0505','fontSize': 20,'textAlign':'center', 'background-color': '#FA0505'})),

]


layout = html.Div([
        dbc.Row(dbc.Col(html.H1("Mapping OMICs to KEGG pathways", style={'textAlign':'center','fontSize': 30}))),
        html.Hr(),
        dbc.Row(dbc.Col(html.H4("Select OMICs and condition data:", style={'fontSize': 16}),
                        width={'size': 3, 'offset': 0}), className="mx-2"),
        dbc.Row([dbc.Col(dcc.Dropdown(id='path_input', value=pathways.Name[1],clearable=False,
                                    persistence=True, persistence_type='memory',
                                    options=[{'label': x, 'value': x} for x in pathways.Name],
                                    style={'fontSize': 14}), 
                                    width={'size': 4, "offset": 0}),
            dbc.Col(dcc.Dropdown(id='PCA_data_options', clearable=False, value='LB60',
                                    persistence=True, persistence_type='memory',
                                    options=[{'label': x, 'value': x} for x in PCA_data_options],
                                    style={'fontSize': 14}), 
                                    width={'size': 2, "offset": 0}),
            dbc.Col(dcc.Dropdown(id='protein_or_transcript', clearable=False, value='Transcriptomics',
                                    options=[{'label': x, 'value': x} for x in t_or_p],
                                    style={'fontSize': 14}), 
                                    width={'size': 2, "offset": 0}),
            
            
            dbc.Col(html.Div([dbc.Button("Color", style={'fontSize': 14},id="hover-target-3", className="me-1", n_clicks=0,),
            dbc.Popover(card_content_colors,target="hover-target-3",placement="right",offset="0,20",body=True,trigger="hover", ),
                ])),
            ], className="mx-2"),
        html.Hr(),
        dbc.Card([dbc.CardBody(id='output-coa',)], color='secondary',),
            
        dbc.Row(
                [
                    dbc.Col(dbc.Card(card_content_gene, color="secondary border-grey", outline=True),),
                     dbc.Col(dbc.Card(card_content_compound, color="secondary border-grey", outline=True)),
                     ],),
])

@app.callback(
    [
        Output('output-coa', 'children'),
             Output('compound-for-kegg', 'children'),
                Output('gene-for-kegg', 'children')],
              [Input('path_input', 'value'),
              Input('PCA_data_options','value'),
              Input('protein_or_transcript', 'value'),])
def show_coa(path_input, PCA_data_options,protein_or_transcript):
    if path_input is not None and PCA_data_options is not None:
        display_path_input = pathways.loc[pathways.Name==path_input,'KEGG_Gene'].item()
        display_path_input = display_path_input[5:]
        if protein_or_transcript == 'Transcriptomics':
            PCA = make_df('all_transcript', countData)
            PCA['LB60']=PCA['60gLNaCL_RPKM']
            PCA['LB20']=PCA['20gLNaCL_RPKM']
            PCA['LB100']=PCA['100gLNaCL_RPKM']
            PCA['HighN']=PCA['Urea3.6gL_RPKM']
            PCA['9hr']=PCA['09Fermentation_RPKM']
            PCA['19hr']=PCA['19Fermentation_RPKM']
            PCA['30hr']=PCA['30Fermentation_RPKM']
            PCA['30 - 9hr']=PCA['30Fermentation_RPKM']-PCA['09Fermentation_RPKM']
            PCA['LB20 - LB60']=-PCA['60gLNaCL_RPKM']+PCA['20gLNaCL_RPKM']
            PCA['LB100 - LB60']=-PCA['60gLNaCL_RPKM']+PCA['100gLNaCL_RPKM']
        elif protein_or_transcript == 'Proteomics':
            PCA = make_df('all_protein', PCA)
            PCA['LB60']=PCA['60gLNaCL_ProteinpgDW']
            PCA['LB20']=PCA['20gLNaCL_ProteinpgDW']
            PCA['LB100']=PCA['100gLNaCL_ProteinpgDW']
            PCA['HighN']=PCA['Urea3.6gL_ProteinpgDW']
            PCA['9hr']=PCA['09Fermentation_ProteinpgDW']
            PCA['19hr']=PCA['19Fermentation_ProteinpgDW']
            PCA['30hr']=PCA['30Fermentation_ProteinpgDW']
            PCA['30 - 9hr'] = PCA['30Fermentation_ProteinpgDW'] - PCA['09Fermentation_ProteinpgDW']
            PCA['LB20 - LB60']=-PCA['60gLNaCL_ProteinpgDW']+PCA['20gLNaCL_ProteinpgDW']
            PCA['LB100 - LB60']=-PCA['60gLNaCL_ProteinpgDW']+PCA['100gLNaCL_ProteinpgDW']
        to_draw_proteomics = PCA[["KEGGSYMBOL",PCA_data_options]]   
        to_draw_proteomics = to_draw_proteomics.iloc[:,0:2]      
        to_draw_proteomics = to_draw_proteomics.rename(columns={PCA_data_options: 'val', 'KEGGSYMBOL': 'KEGGSYMBOL'})

        pathway = KGML_parser.read(kegg_get(display_path_input, "kgml"))
        pathstring_compound, pathstring_gene = its_a_kegg_pathway(display_path_input)

        bin_labels_8 = ['#060EFB', '#4449F0', '#B2AAF1', '#c4fdf0', '#fbfb8d', '#FDA88C', '#FC5749', '#FA0505']
        to_draw_proteomics['eight_bins'] = pd.cut(to_draw_proteomics.val, bins=8, labels=bin_labels_8)
        to_draw_proteomics['eight_bins_cts'] = pd.cut(to_draw_proteomics.val, bins=8)
        to_draw_proteomics = to_draw_proteomics[to_draw_proteomics['eight_bins'].str.contains('#')]
        for element in pathway.orthologs:
            for graphic in element.graphics:
                graphic.fgcolor = '#ffffff'
                graphic.bgcolor = '#ffffff'
        for element in pathway.orthologs:
            for graphic in element.graphics:
                for ko in to_draw_proteomics.KEGGSYMBOL:
                    if ko in element.name:
                        graphic.fgcolor = to_draw_proteomics.loc[to_draw_proteomics.KEGGSYMBOL == ko,'eight_bins'].values[0]
                        graphic.bgcolor = to_draw_proteomics.loc[to_draw_proteomics.KEGGSYMBOL == ko,'eight_bins'].values[0]   
    #f = df['eight_bins_cts'].value_counts(sort=False).reset_index()
    #return df['eight_bins_cts'].value_counts(sort=False).reset_index()
       
        canvas = KGMLCanvas(pathway, import_imagemap=True)
        buffered = io.BytesIO()
        canvas.draw(buffered)
        new = base64.b64encode(buffered.getvalue())
        encoded = bytes("data:application/pdf;base64,", encoding='utf-8') + new


        return html.Div([html.Img(src=encoded.decode('utf-8')),html.Hr(),]), pathstring_compound, pathstring_gene