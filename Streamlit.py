#!/usr/bin/env python
# coding: utf-8

# In[1]:


import streamlit as st
st.set_page_config(layout="wide")
import Bio
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import base64
from IPython.display import Image, HTML, display
import pandas as pd
import numpy as np
import os
import altair as alt
from sklearn import decomposition
st.set_option('deprecation.showPyplotGlobalUse', False)
from Bio.SeqUtils.IsoelectricPoint import IsoelectricPoint as IP
matplotlib.pyplot.switch_backend('Agg') 

# cobra / escher
from escher import Builder
import cobra
from cobra import Metabolite, Reaction, summary
from cobra.io import read_sbml_model, write_sbml_model, save_json_model, load_json_model
from cobra.flux_analysis import (
    single_gene_deletion, single_reaction_deletion, double_gene_deletion,
    double_reaction_deletion, production_envelope)

def to_df(result):  #make nice df
    return pd.read_table(io.StringIO(result),names=['KEGG Gene','Name'])
def PDF(filename): # not currently used
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

def head(text, lines=10): # A bit of helper code to shorten long text
    return text.split("\n")[:lines]

def draw_kegg_map(map_id, outdir, df): #kegg map by binning into colors
    pathway = KGML_parser.read(kegg_get(map_id, "kgml"))
    #yellow to red #bin_labels_8 = ['#fbfb04', '#ffe100', '#ffc600', '#ffab00', '#ff8d00', '#ff6e00', '#ff4900', '#fc020c']
    bin_labels_8 = ['#060EFB', '#4449F0', '#B2AAF1', '#c4fdf0', '#fbfb8d', '#FDA88C', '#FC5749', '#FA0505']
    df['eight_bins'] = pd.cut(df.val, bins=8, labels=bin_labels_8)
    df['eight_bins_cts'] = pd.cut(df.val, bins=8)
    df = df[df['eight_bins'].str.contains('#')]
    for element in pathway.orthologs:
        for graphic in element.graphics:
            graphic.fgcolor = '#ffffff'
            graphic.bgcolor = '#ffffff'
    for element in pathway.orthologs:
        for graphic in element.graphics:
            for ko in df.KEGGSYMBOL:
                if ko in element.name:
                    graphic.fgcolor = df.loc[df.KEGGSYMBOL == ko,'eight_bins'].values[0]
                    graphic.bgcolor = df.loc[df.KEGGSYMBOL == ko,'eight_bins'].values[0]   
    
    canvas = KGMLCanvas(pathway, import_imagemap=True)
    img_filename = "%s.pdf" % map_id
    canvas.draw(os.path.join(outdir, img_filename))
    f = df['eight_bins_cts'].value_counts(sort=False).reset_index()
    return df['eight_bins_cts'].value_counts(sort=False).reset_index()

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

def make_df(choice):
    if choice is 'all_data': 
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLProteinpgDW', '20gLProteinpgDW','NaCl100g/LProteinpgDW','Urea3.6g/LProteinpgDW','Fermentation9hProteinpgDW', 'Fermentation19hProteinpgDW','Fermentation30hProteinpgDW', 'NaCl60gLRPKM', 'NaCl20gLRPKM','NaCl100g/LRPKM','Urea3.6g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM', 'Fermentation30hRPKM']].copy()
    elif choice is 'fermentation_protein':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','Fermentation9hProteinpgDW', 'Fermentation19hProteinpgDW', 'Fermentation30hProteinpgDW']].copy()
    elif choice is 'fermentation_transcript':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','Fermentation9hRPKM','Fermentation19hRPKM', 'Fermentation30hRPKM']].copy()
    elif choice is 'fermentation_all':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','Fermentation9hRPKM','Fermentation19hRPKM', 'Fermentation30hRPKM','Fermentation9hProteinpgDW', 'Fermentation19hProteinpgDW','Fermentation30hProteinpgDW']].copy()
    elif choice is 'conditions_protein':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLProteinpgDW', '20gLProteinpgDW','NaCl100g/LProteinpgDW', 'Urea3.6g/LProteinpgDW']].copy()
    elif choice is 'conditions_transcript':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','NaCl60gLRPKM', 'NaCl20gLRPKM','NaCl100g/LRPKM']].copy()
    elif choice is 'all_transcript':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','NaCl60gLRPKM', 'NaCl20gLRPKM','NaCl100g/LRPKM','Urea3.6g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM', 'Fermentation30hRPKM']].copy()
    elif choice is 'all_protein':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLProteinpgDW', '20gLProteinpgDW','NaCl100g/LProteinpgDW','Urea3.6g/LProteinpgDW','Fermentation9hProteinpgDW', 'Fermentation19hProteinpgDW','Fermentation30hProteinpgDW']].copy()        
    elif choice is 'conditions_all':
        PCA = countData[['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','NaCl60gLRPKM', 'NaCl20gLRPKM','NaCl100g/LRPKM', 'Urea3.6g/LRPKM','60gLProteinpgDW', '20gLProteinpgDW','NaCl100g/LProteinpgDW', 'Urea3.6g/LProteinpgDW']].copy()

    for column in PCA.columns.drop(['Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL']):
        PCA[column] = standardize_2(PCA[column])
    return PCA
# Cobrapy functions
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
def data_frame_output_flux(model,name):
    a=model.optimize()
    b=model.summary()
    fluxes_df=pd.DataFrame(b.uptake_flux)
    c=pd.DataFrame(b.secretion_flux)
    fluxes_df=pd.concat([fluxes_df,c])
    fluxes_df=fluxes_df.rename(columns={'flux':name})
    return fluxes_df.drop(['reaction','metabolite'],axis=1)
def make_chart_varying_media(model,glcs,kind,media,upper=True,both=True):
    with model:
        rgs = []
        for glc in glcs:
            if both:
                model.reactions.get_by_id(media).lower_bound = glc
                model.reactions.get_by_id(media).upper_bound = glc
            elif upper:
                model.reactions.get_by_id(media).upper_bound = glc
            else:
                model.reactions.get_by_id(media).lower_bound = glc
            if kind == 'Growth':
                s3 = model.optimize()
            elif kind == 'Ectoine':
                s3 = test_ectoine_production(model)
            elif kind == 'PHA':
                s3 = test_PHA_production(model)
            if s3.objective_value is None:
                rgs.append(0)
            else:
                rgs.append(s3.objective_value)
        return rgs
def testing_amino_acids(lb60,lb20,lb100,name,media,lb=0,ub=1000,a=50,kind='Growth',
                        upper=True,both=True):
    plt.clf()
    glcs = np.arange(lb,ub,a)
    df = pd.DataFrame(columns={name,'LB60','LB20','LB100'})
    df[name] = glcs
    df.LB60 = make_chart_varying_media(lb60,glcs,kind,media,upper,both)
    df.LB100 = make_chart_varying_media(lb100,glcs,kind,media,upper,both)
    df.LB20 = make_chart_varying_media(lb20,glcs,kind,media,upper,both)
    
    df = df.melt(name, var_name='media',value_name=kind)
    return sns.lineplot(x=name,y=kind,hue="media",data=df)

amino_acids = ['CYS_c','L-ASPARTATE_c','GLT_c','PHE_c','GLY_c','HIS_c','ILE_c','LYS_c','LEU_c','MET_c','ASN_c','PRO_c','GLN_c','ARG_c','SER_c','THR_c','VAL_c','TRP_c','L-ALPHA-ALANINE_c','TYR_c','ADENINE_c']
def test_aa_production(model):
    rgs = pd.DataFrame(columns = ['Output'])
    i=0
    for met_id in amino_acids:
        aa_sink = cobra.Reaction('aa_sink')
        aa_sink.add_metabolites({model.metabolites.get_by_id(met_id):-1})
        with model:
            model.add_reactions([aa_sink]) 
            set_yeast_extraction(model,lb=0,ub=0)
            model.objective = 'aa_sink'
            model.objective_direction = 'max'
            s1 = model.optimize()
    
        rgs.at[i,'Output'] = s1.objective_value
        i+=1
    return rgs
carbon_sources = ['Exchange_Glucopyranose','FRUtex','Exchange_L-LACTATE','Exchange_ACET','Exchange_CIT','GLYCtex','Exchange_Trehalose','BUTtex','EX_SUC_e','ETOHtex']
def reset_carbon(model,input_for_conditions,yeast_adjust):
    for c in carbon_sources:
        model.reactions.get_by_id(c).upper_bound = 0
    model.reactions.get_by_id(input_for_conditions).upper_bound = 0
    if yeast_adjust == 'Minimal':
        set_yeast_extraction(model,0,0)
    elif yeast_adjust == 'LB':
        set_yeast_extraction(model,10,-10)
    elif yeast_adjust == 'LB+':
        set_yeast_extraction(model,100,-100)
    
#import the data
read_and_cache_csv = st.cache(pd.read_csv)
countData = read_and_cache_csv('https://raw.githubusercontent.com/helloftroy/Omics_TD1.0/main/R_Landing.csv')
scatterData = read_and_cache_csv('https://raw.githubusercontent.com/helloftroy/Omics_TD1.0/main/Scatter_New.csv')
countData = pd.DataFrame(countData)
#make shorter pathway list
#pathlist_2 = ['Signal transduction','Transport and catabolism','Cell growth and death','Cellular community - eukaryotes','Cell motility','Membrane transport','Replication and repair','Folding, sorting and degradation','Translation','Transcription','Xenobiotics biodegradation and metabolism','Biosynthesis of other secondary metabolites','Metabolism of terpenoids and polyketides','Metabolism of cofactors and vitamins','Glycan biosynthesis and metabolism','Oxidative phosphorylation','Photosynthesis','Photosynthesis - antenna proteins','Carbon fixation in photosynthetic organisms','Carbon fixation pathways in prokaryotes','Glyoxylate and dicarboxylate metabolism','Propanoate metabolism','Butanoate metabolism','C-Branched dibasic acid metabolism','Metabolism of other amino acids','Amino acid metabolism','Quorum sensing','Nucleotide metabolism','Lipid metabolism','Fatty acid degradation','Inositol phosphate metabolism','Methane metabolism','Nitrogen metabolism','Sulfur metabolism','Glycolysis / Gluconeogenesis','Citrate cycle (TCA cycle)','Pentose phosphate pathway','Pentose and glucuronate interconversions','Fructose and mannose metabolism','Galactose metabolism','Ascorbate and aldarate metabolism','Starch and sucrose metabolism','Amino sugar and nucleotide sugar metabolism','Pyruvate metabolism']
pathlist_2 = ['Ribosom', 'Carbon fixation pathways in prokaryotes', 'ABC transporter', 'Glyoxylate and dicarboxylate metabolism', 'Arginine biosynthesis', 'Butanoate metabolism', 'Quorum sensing',  'Fatty acid metabolism',  'Citrate cycle (TCA cycle)', 'Pyrimidine metabolism', '2-Oxocarboxylic acid metabolism', 'Synthesis and degradation of ketone bodies', 'Biosynthesis of amino acids','Glycine, serine and threonine metabolism', 'Protein export', 'Pentose phosphate pathway', 'Arginine and proline metabolism',  'Pyruvate metabolism', 'Oxidative phosphorylation','Phosphotransferase system (PTS)',  'DNA replication',   'Glutathione metabolism',    'Tyrosine metabolism', 'Lysine degradation','Homologous recombination', 'Folate biosynthesis', 'Galactose metabolism','Phenylalanine metabolism', 'RNA polymerase', 'Biofilm formation - Escherichia coli', 'Flagellar assembly',  'RNA degradation', 'C5-Branched dibasic acid metabolism', 'Bacterial chemotaxis', 'RNA transport', 'Vitamin B6 metabolism','Purine metabolism', 'Propanoate metabolism', 'Lipoic acid metabolism', 'Glycolysis / Gluconeogenesis','Selenocompound metabolism', 'Aminoacyl-tRNA biosynthesis','Amino sugar and nucleotide sugar metabolism','Lipopolysaccharide biosynthesis', 'Base excision repair','Nicotinate and nicotinamide metabolism', 'Ubiquinone and other terpenoid-quinone biosynthesis','Alanine, aspartate and glutamate metabolism', 'Glycerophospholipid metabolism', 'Thiamine metabolism',  'Peptidoglycan biosynthesis', 'Mismatch repair','Phenylalanine, tyrosine and tryptophan biosynthesis', 'Sulfur metabolism', 'Cysteine and methionine metabolism', 'Ascorbate and aldarate metabolism', 'Methane metabolism', 'Biosynthesis of unsaturated fatty acids','Terpenoid backbone biosynthesis', 'Riboflavin metabolism', 'Sulfur relay system',  'D-Glutamine and D-glutamate metabolism','Nucleotide excision repair', 'Thermogenesis', 'Pyrimidine metabolism', 'One carbon pool by folate', 'Nitrogen metabolism', 'Glycerolipid metabolism', 'Inositol phosphate metabolism','Glucosinolate biosynthesis',       'Fructose and mannose metabolism', 'Valine, leucine and isoleucine biosynthesis',  'Valine, leucine and isoleucine degradation', 'Fatty acid degradation']

countData['KEGG_Pathway'] = 'Missing'
i=0
for key in countData.KEGGPATH:
    for kegginput in pathlist_2:
        if kegginput in key:
            countData.loc[i,'KEGG_Pathway'] = kegginput
    i=i+1
lb60_model = read_sbml_model('/Users/hellpark/Desktop/Bioinformatics September 2021/model_objects/Chalmers_Models/lb60_model.sbml')
lb20_model = read_sbml_model('/Users/hellpark/Desktop/Bioinformatics September 2021/model_objects/Chalmers_Models/lb20_model.sbml')
lb100_model = read_sbml_model('/Users/hellpark/Desktop/Bioinformatics September 2021/model_objects/Chalmers_Models/lb100_model.sbml')
    
###
###
###
###
### Begin Main Code
###
###
###
choice = st.sidebar.radio("Chart: ",('PCA-George','KEGG', 'Cobrapy Cell Model'))
outdir="/Users/hellpark/Desktop/" 

if choice == 'PCA-George':
    
    choice_3 = st.sidebar.radio("Remove genes w/o KEGG pathways?", ('No','Yes'))
    st.write("# PCA (Principle Component Analysis)")

    # PCA 1
    PCA = make_df('all_data')
    PCA = PCA.dropna()
    if choice_3 == 'Yes':
        PCA = PCA[PCA.KEGG_Pathway != 'Missing']
    df_genes = PCA.drop(labels={'Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL'}, axis=1)
    to_plot_and_analyze, loadings, expl_var = PCA_script(df_genes, transpose=True, loadings = True, pcs = 3)
    to_plot_and_analyze = pd.DataFrame(to_plot_and_analyze, columns = ['PC1', 'PC2', 'PC3'])
    to_plot_and_analyze['Condition'] = ['LB60','LB20','LB100','High N','9Hr','19Hr','30Hr','LB60','LB20','LB100','High N','9Hr','19Hr','30Hr']
    to_plot_and_analyze['Labels'] = ['LB60 Protein','LB20 Protein','LB100 Protein','High N Protein','9Hr Protein','19Hr Protein','30Hr Protein','LB60 Transcript','LB20 Transcript','LB100 Transcript','High N Transcript','9Hr Transcript','19Hr Transcript','30Hr Transcript']
    PCA['loadings PC1'] = loadings[:,0]
    PCA['loadings PC2'] = loadings[:,1]
    PCA['loadings PC3'] = loadings[:,2]
    
    st.markdown("PCA is an unsupervised machine-learning method that aims to reduce dimensionality. PCA can take in large data sets, and through linear transformation, simplify while minimizing information loss. Therefore, we can learn which indicators drive the most variance between different conditions. [1]")
    st.markdown("Figure 1 shows the initial PCA for 7 conditions, Transcriptomics and Proteomics (~600 non-zero proteins).")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown(" - Current standardization method: Take Log2 of all non-NaN rows. Take (x-mean(column)) / stdev(column).")
        st.markdown(" - This is standard practice to make sure the different data sources are all normalized to be same scale. ")
        st.markdown('#')
        to_plot_and_analyze['abs(PC1)'] = np.absolute(to_plot_and_analyze['PC2'])
        chart = alt.Chart(to_plot_and_analyze).mark_circle(size=200).encode(x='PC1',y='PC2',
                                            color='Condition',tooltip='Labels')
        st.markdown("    #### Fig 1) PCA for all 7 conditions")
        st.altair_chart(chart)    
    with col2:
        st.dataframe(pd.DataFrame(expl_var, index=['Explained Variance PC1','Explained Variance PC2','Explained Variance PC3'], columns=['%']))
        st.markdown('#')
        chart_2 = alt.Chart(to_plot_and_analyze).mark_circle(size=200).encode(x='abs(PC1)',y='PC2',
                                            color='Condition',tooltip='Labels')
        #st.markdown("""<div style="text-align: center"> Fig 2) PCA with absolute value of PC1 </div>""",unsafe_allow_html=True)
        st.markdown("    #### Fig 2) PCA with absolute value of PC1")
        st.altair_chart(chart_2)
    st.markdown("## Thoughts:")
    st.markdown("1) It is super that have 75% variance captured first two PCs, means that 72% variance is not lost.")
    st.markdown("2) PC1's output is clear but surprising: the data type (Transcript/Proteomic) explain 55% of the variation. This means there are clear differences in expression between what is transcribed and what is translated across **ALL the conditions**. In other words, DNA and Protein are strongly anticorrelated in TD1.0.") 
    st.markdown("3) PC2's explains 20% of the variation. Upon further investigation, PC2 appears to separate based on experimental condition. This becomes obvious when taking the absolute value of PC1 (Figure 2.) All points line up, except a slight discrepancy between High Nitrogen and the 9Hr fermentation conditions. The consistancy is an indication our PCA approach is validated to be working properly.")
    st.markdown("#### Three questions arise: **What drives the differences between conditions?**, **What drives the differences between transcription, and translation?** and **Why do the DNA and protein variables pair up so nicely?** In the next section, we will try to enhance those differences.")
    
                
    st.markdown("## Corrections for PC1")
    st.markdown("Since PC1 is a clear separation between DNA/Protein, our next step is to correct for this pattern. We rotate the data matrix prior to PCA, via the following transformation; **A1** is magnitude based, so will enhance **differences between conditions** and **A2** is subtraction based, and should amplify the **differences between DNA and Protein** for further exploration.")
    st.markdown("- **A1**: (DNA_i + Protein_i) and **A2**: (DNA_i - Protein_i), for i = 1,2,3,4,5,6,7. ")    
    #transformed PCA
    PCA = make_df('all_data')
    PCA = PCA.dropna()
    if choice_3 == 'Yes':
        PCA = PCA[PCA.KEGG_Pathway != 'Missing']
    
    PCA['LB60+']=PCA['60gLProteinpgDW']+PCA['NaCl60gLRPKM']
    PCA['LB20+']=PCA['20gLProteinpgDW']+PCA['NaCl20gLRPKM']
    PCA['LB100+']=PCA['NaCl100g/LProteinpgDW']+PCA['NaCl100g/LRPKM']
    PCA['HighN+']=PCA['Urea3.6g/LProteinpgDW']+PCA['Urea3.6g/LRPKM']
    PCA['9hr+']=PCA['Fermentation9hProteinpgDW']+PCA['Fermentation9hRPKM']
    PCA['19hr+']=PCA['Fermentation19hProteinpgDW']+PCA['Fermentation19hRPKM']
    PCA['30hr+']=PCA['Fermentation30hProteinpgDW']+PCA['Fermentation30hRPKM']
    PCA['60-']=PCA['60gLProteinpgDW']-PCA['NaCl60gLRPKM']
    PCA['20-']=PCA['20gLProteinpgDW']-PCA['NaCl20gLRPKM']
    PCA['100-']=PCA['NaCl100g/LProteinpgDW']-PCA['NaCl100g/LRPKM']
    PCA['N-']=PCA['Urea3.6g/LProteinpgDW']-PCA['Urea3.6g/LRPKM']
    PCA['9-']=PCA['Fermentation9hProteinpgDW']-PCA['Fermentation9hRPKM']
    PCA['19-']=PCA['Fermentation19hProteinpgDW']-PCA['Fermentation19hRPKM']
    PCA['30-']=PCA['Fermentation30hProteinpgDW']-PCA['Fermentation30hRPKM']

    df_genes = PCA.drop(labels={'Fasta','Name','ECNumber','KEGG_Pathway','KEGGSYMBOL','60gLProteinpgDW','20gLProteinpgDW', 'NaCl100g/LProteinpgDW', 'Urea3.6g/LProteinpgDW','Fermentation9hProteinpgDW', 'Fermentation19hProteinpgDW','Fermentation30hProteinpgDW', 'NaCl60gLRPKM', 'NaCl20gLRPKM','NaCl100g/LRPKM', 'Urea3.6g/LRPKM', 'Fermentation9hRPKM','Fermentation19hRPKM', 'Fermentation30hRPKM'}, axis=1)
    to_plot_and_analyze, loadings, expl_var = PCA_script(df_genes, transpose=True, loadings = True, pcs = 3)
    to_plot_and_analyze = pd.DataFrame(to_plot_and_analyze, columns = ['PC1', 'PC2', 'PC3'])
    to_plot_and_analyze['Condition'] = ['LB60','LB20','LB100','N','9Hr','19Hr','30Hr','LB60','LB20','LB100','N','9Hr','19Hr','30Hr']
    to_plot_and_analyze['Labels'] = ['LB60 added','LB20 added','LB100 added','N added','9Hr added','19Hr added','30Hr added','LB60 subtracted','LB20 subtracted','LB100 subtracted','N subtracted','9Hr subtracted','19Hr subtracted','30Hr subtracted']
    PCA['loadings PC1'] = loadings[:,0]
    PCA['loadings PC2'] = loadings[:,1]
    PCA['loadings PC3'] = loadings[:,2]
    
    col1, col2 = st.columns(2)
    with col1:
        chart_3 = alt.Chart(to_plot_and_analyze).mark_circle(size=200).encode(x='PC1',y='PC2',
                                            color='Condition',tooltip='Labels')
        st.dataframe(pd.DataFrame(expl_var, index=['Explained Variance PC1','Explained Variance PC2','Explained Variance PC3'], columns=['%']))
        st.markdown("    #### Fig 3) Rotated matrix PCA")
        st.altair_chart(chart_3)


    with col2:
        st.markdown("#### Transforming matrix to reduce DNA and Protein effects:")
        st.markdown("- An initial check at our technique indicates that the rotated matrix still shows segregation along PC1, which explains 80% of data variation. ")
        st.markdown("- In PC2, the **A2**, (DNA - Protein) case clusters together. This indicates that there is little variation between conditions in terms of DNA and Protein expression, in other words **differences between DNA and Protein do not drive the variation between each experimental condition.** Variation between DNA/Protein exist in all conditions.")
        st.markdown("- Meanwhile, in **A1**, (DNA + Protein) variation between conditions is enhanced. Specifically, LB20 is extremely ++ in PC2, LB100 is also +, and then the other conditions are all -, with the 30-hour fermentation most different from LB20. ")
    st.markdown("### So which genes cause most variation in the data?")
    st.markdown("It is logical to now separate these two observations to understand what is driving both. We can do this by performing PCA on ONLY **A1** and **A2** and perform a pathway and gene cluster analysis")          
    st.markdown("- We can look at the loadings for each gene. The loadings are like weights: a high loading indicates that that gene drives the variation, while a loading near zero means it drives little variation.")
    st.markdown("- Recall anything very positive describes LB20, LB100. Negative describes anything closer to later fermentation samples.")
    col1, col2 = st.columns(2)
    with col1:
        largest = PCA.nlargest(5, 'loadings PC2')
        st.markdown("""<style>tbody th {display:none}.blank {display:none}</style>""", unsafe_allow_html=True)
        st.table(largest.iloc[:,[3,1,34]])
    with col2:
        smallest = PCA.nsmallest(5, 'loadings PC2')
        st.markdown("""<style>tbody th {display:none}.blank {display:none}</style>""", unsafe_allow_html=True)
        st.table(smallest.iloc[:,[3,1,34]])
    
    st.markdown("# Observations for most-influential genes")
    st.markdown("1) In the top 10 with largest negative loading (upregulated in LB20) 5 are related to flagellar assembly and cell motility, while the other 5 are related to cell maintainence: RNA polymerase, tRNA synthase, or ribosomal proteins.")
    st.markdown("- **FliH, FliF, and FliM** are all tied to **flagellar assembly and chemotaxis.**")
    st.markdown("- **The PTS (phototransferase system)** membrane protein's main function is to translocate sugars across the membrane whilst phosphorylation. In L. monocytogenes, PTS fructose-specific protein was tied to improved salt stress resistance [6]. Specifically, cold-stressed cells overexpressed this protein, and were then more protected when exposed to salt. **Are fructose, mannose and other sugars are important osmoprotectants?** PTS enzyme IIABC can convert glycoproteins to be a compatible solute, and PTS has previously been tied to salt concentration. **ATP-dependent RNA helicase** may be indirectly ties to **biofilm formation.** E.coli mutants have reduced biofilm, and increased cell dispersion. Moreover, KO's proteome indicated a downregulation in flagellin paired with slower observed motility [7]. The pseudoridine-55 synthase gene (truB in E.coli), has also been tied to **outer membrane biogenesis.** Specifically, truB mutants also had abnormal levels of outer membrane proteins OmpA and OmpX, causing competitive disadvantage but not cell death.")    
    st.markdown("2) In the top 10 largest positive loading, (corresponding to downregulation in LB20) many ABC transporters are listed, including transporters tied to quorum sensing. the short-chain oxidoreductase/dehydrogenase catalyzes glucose & NAD+ -> NADH, and has been reported to be upregulated in B. pseudomallei until salt stress [8].")
    st.markdown("#### Try switching off missing KEGG pathways - the PCA results look similar but the proteins analyzed are all annotated. Upon first look, the most ++ genes are related to cell motility, and PTS system. The most -- are ABC transporters, or involved in lipid metabolism.")
    
    st.markdown("## Cluster Heat Map")
    st.markdown("Our goal is to detect broad patterns of gene expression for differences between conditions. Genes are at right. At the bottom 'T' are transcriptomics, 'P' are proteomics. Red=high, blue=low expression. **Again, can try turning off missing KEGG pathways at left for more meaningful results.**")
    heatmap = PCA.nlargest(40, 'loadings PC2')
    heatmap = heatmap.append(PCA.nsmallest(40, 'loadings PC2'))
    label_heatmap = ['LB60 P','LB20 P','LB100 P','High N P','9Hr P','19Hr P','30Hr P','LB60 T','LB20 T','LB100 T','High N T','9Hr T','19Hr T','30Hr T']
    ylabels = heatmap['Name'].str[:28]
    network_pal = sns.husl_palette(14, s=.45)
    network_lut = dict(zip(map(str, label_heatmap), network_pal))
    sns.set(font_scale=0.6)
    heatmap_plot = sns.clustermap(heatmap[['60gLProteinpgDW', '20gLProteinpgDW','NaCl100g/LProteinpgDW','Urea3.6g/LProteinpgDW','Fermentation9hProteinpgDW', 'Fermentation19hProteinpgDW','Fermentation30hProteinpgDW', 'NaCl60gLRPKM','NaCl20gLRPKM','NaCl100g/LRPKM','Urea3.6g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM', 'Fermentation30hRPKM']],xticklabels=label_heatmap, yticklabels=ylabels,method='ward',standard_scale=1,cmap="coolwarm",figsize=(7,9))
    st.pyplot(heatmap_plot)
    st.markdown("The heatmap above grabs the top genes with most positive and most negative loadings in PC2. A few observations stand out:")
    st.markdown("- LSU/SSU ribosomal proteins tend to be highly expressed in all protein conditions, especially at 9hr fermentation. Noticeably, ribosomal proteins are also active in 20/100 g/L salt transcripts but not in the normal salt conditions.")
    st.markdown("- There is higher transporter activity for normal salt than for high/low salt in proteomics data. LB20 stands out with much lower transporter activity. ")

    st.markdown("## Pathway Enrichment Analysis")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("The heatmap showed no obvious patterns, next perform enrichment analysis of all genes, using PCA loading, and binned using KEGG IDs to KEGG pathways. Below toggle between average of loadinings, and sum of loadings. **Which KEGG pathways are the most different between conditions?**")
    with col2:
        choice_4 = st.radio("Loadings calculation", ('Average','Sum'))
    if choice_4 == 'Average':
        PCA_KEGG_grouped = PCA.groupby('KEGG_Pathway').mean().reset_index()
    elif choice_4 == 'Sum':
        PCA_KEGG_grouped = PCA.groupby('KEGG_Pathway').sum().reset_index()
    chart_4 = alt.Chart(PCA_KEGG_grouped).mark_bar().encode(y=alt.Y('KEGG_Pathway', sort='-x'),x='loadings PC2')
    st.altair_chart(chart_4,use_container_width=True)
    
    st.markdown("## Initial observations for pathway analysis: Upregulated Pathways")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### Flagella and cell motility upregulation")
        st.markdown(" The next most differentiated pathways are in flagella assembly and cell motility. **Flagellar are large and complex molecules, why would sick cells devote energy to them?**")
        st.markdown("- TD1.0 might be using flagella to propell forward, into better conditions. Especially at low osmotic pressure, swimming motility increases in some bacteria like E. albertii [9]. The flagella gene (RflP) in E.coli has been shown to responds to cell envelope stress, alterations of outer membrane integrity. Flagella can help cells survive stress; E.coli cells survives higher concentrations of gentamicin when expressing flagella. When add combatible solutes such as glutamic acid or proline, the swimming cells decreased 35-fold, suggesting osmotic pressue can be a trigger of flagella in some organisms.")
        st.markdown("- The inherent trade-off between growth and protein production suggests that resources for protein production may be liberated following growth-inhibiting bioprocess treatments. SO won't make enough ectoine/PHA if everything allocated to flagella")    
    with col2:
        st.markdown("#### Ribosome upregulation")
        st.markdown(" Previously, [22] found that ribosomal proteins were increased at high salt. We see this too (LB100>LB60), but we also see high ribosomal proteins, the highest, at low salt. This indicates the cells upregulate protein biosynthesis when outside typical salt range (check ko03010 in KEGG tab.) In theory, moderate halophiles should have high enough concentrations of organic compatible solutes, and their proteins should not denature or require replenishment. However, moderate halophiles still have slightly acidic proteomes; at high salt, if protein aggregation threatens halobacteria can exhibit excess acidic/deficite basic amino acids on their protein surfaces. Excess negative charge have been shown to be localized to the surface of modelled proteins [20]. **What is the ribosome translating? Where are cellular resources being allocated?**")
   
    PCA2 = make_df('conditions_transcript')
    PCA2 = PCA2.dropna()

    PCA2['LB20 - LB60']=-PCA2['NaCl60gLRPKM']+PCA2['NaCl20gLRPKM']
    to_draw_proteomics = PCA2[["KEGGSYMBOL",'LB20 - LB60']]
    to_draw_proteomics = to_draw_proteomics.rename(columns={'LB20 - LB60': 'val', 'KEGGSYMBOL': 'KEGGSYMBOL'})
    img_filename = "%s.pdf" % 'ko02040'
    with open(os.path.join(outdir, img_filename),"rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="900" height="700" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)
        
    st.markdown("## Initial observations for pathway analysis: Downregulated pathways")
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("#### Butanoate Metabolism")
        st.markdown("Butanoate metabolism is downregulated in LB20, but no other condition. Notably, this pathway contains the phaCAB cassette. It seems PHA synthesis is downregulated in LB20. phaA, phaC, phaB are all decreased along with adjacent enzymes. Notably, LB100 has similar expression to LB60.")
        st.markdown("Perhaps, this indicates **PHA is used in osmotic stress.** PHA are carbon/energy sources (this is well understood) and can be degraded in unbalanced/non-optimal growth conditions. Wang et al. found that in H mediterranei salt higher than optimal growth conditions led to higher PHA, concluding that lower salt prevents PHA synthesis. [2]")                    
    with col2:
        st.markdown("#### Transporter activity")
        st.markdown("When use sum of pathways, ABC transport is most downregulated in LB20 (try ko02010 in KEGG tab.) It is well known osmolytes can be used by halotolerant species to combat osmotic pressure. However, each species appears to have its own strategy. For example, some species import proline, others produce ectoine (TD1.0 is known for this.) **We can hypothesis TD1.0's strategies to combat salt by observing which ABC Transporters shift.**") 
        st.markdown("- Comparing LB20 to LB60: Urea and amino acid transporters are extremely downregulated. Glucose, mannose, glycerol, trehalose are all downregulated, but **fructose is heavily upregulated.** Finally, Sec-SRP, an ATPase (producing ATP and transporting in amino acids) is upregulated. Since amino acid transport is decreased, perhaps LB20 halomonas selects for this transporter, overall reducing other amino acid transport (ABC's consume ATP,) reducing overall import but saving energy.")
        st.markdown("- Comparing LB100 to LB60: Glycine/betaine, proline, maltose, trehalose are all upregulated. Glucose is downregulated, but fructose is upregulated. Most amino acids are upregulated, as is iron transport (known activator of ectoine cassette.)")
        st.markdown("This suggests a shift towards **proline, glycine/betaine, and perhaps trehalose** for osmotic protection. The **preference for fructose is interesting**, fructose and mannose metabolism in general are upregulated. ")
   
    value_cnt = draw_kegg_map('ko00650', outdir, to_draw_proteomics)
    img_filename = "%s.pdf" % 'ko00650'
    with open(os.path.join(outdir, img_filename),"rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="900" height="700" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)                            
    
    st.markdown("#### Other Pathways to explore in KEGG tab:")
    st.markdown("- **ko00260:** As expected, there is an overall downregulated in ectoine biosynthesis in LB20 vs LB60, and an overall upregulation in LB100 vs LB60. Since ectoine is a well-known osmoprotectant, it is clear the higher salt necessetates higher ectoine production.")
    st.markdown("- **ko02024:** quorum sensing (quite upregulated in LB100) and ko02030: Bacterial chemotaxis")
    st.markdown("- **ko00330:** Purine and argenine metabolism. Upregulated towards Proline in LB100")
    st.markdown("- **ko00051:** Fructose and mannose metabolism. Upregulated, especially in LB20.")
    st.markdown("- ATP-synthase genes are all upregulated in both LB100 and LB20, indicating vigorous growth and ATP usage is needed for these conditions. Out of 64 annotated LSU/SSU ribosomal proteins, all 64 were upregulated in LB20. ")



   
    
elif choice == 'KEGG':
    #transformed PCA again -> this time can select between proteomics and transcriptomics
    select = st.sidebar.radio("Protein/Transcript", ('Transcript','Protein'))
    if select == 'Transcript':
        PCA = make_df('all_transcript')
        PCA = PCA.dropna()
        PCA['LB60']=PCA['NaCl60gLRPKM']
        PCA['LB20']=PCA['NaCl20gLRPKM']
        PCA['LB100']=PCA['NaCl100g/LRPKM']
        PCA['HighN']=PCA['Urea3.6g/LRPKM']
        PCA['9hr']=PCA['Fermentation9hRPKM']
        PCA['19hr']=PCA['Fermentation19hRPKM']
        PCA['30hr']=PCA['Fermentation30hRPKM']
        PCA['LB20 - LB60']=-PCA['NaCl60gLRPKM']+PCA['NaCl20gLRPKM']
        PCA['LB100 - LB60']=-PCA['NaCl60gLRPKM']+PCA['NaCl100g/LRPKM']
    elif select == 'Protein':
        PCA = make_df('all_protein')
        PCA = PCA.dropna()
        PCA['LB60']=PCA['60gLProteinpgDW']
        PCA['LB20']=PCA['20gLProteinpgDW']
        PCA['LB100']=PCA['NaCl100g/LProteinpgDW']
        PCA['HighN']=PCA['Urea3.6g/LProteinpgDW']
        PCA['9hr']=PCA['Fermentation9hProteinpgDW']
        PCA['19hr']=PCA['Fermentation19hProteinpgDW']
        PCA['30hr']=PCA['Fermentation30hProteinpgDW']
        PCA['LB20 - LB60']=-PCA['60gLProteinpgDW']+PCA['20gLProteinpgDW']
        PCA['LB100 - LB60']=-PCA['60gLProteinpgDW']+PCA['NaCl100g/LProteinpgDW']
        
    conditions = ['LB60','LB20','LB100','HighN','9hr','19hr','30hr', 'LB20 - LB60', 'LB100 - LB60']
    input_for_conditions = st.sidebar.selectbox('Select condition to map:', conditions)
     
    st.write("### List of Halo KEGG pathways")
    st.write(to_df(kegg_list('pathway', 'ko').read()))
    outdir = "/Users/hellpark/Desktop/"
    pathinput = st.text_input("Pathway KEGG (https://www.genome.jp/kegg/pathway.html)",value='ko00010',help='type pathway like, ko00010 after the path:')

    to_draw_proteomics = PCA[["KEGGSYMBOL",input_for_conditions]]
    to_draw_proteomics = to_draw_proteomics.rename(columns={input_for_conditions: 'val', 'KEGGSYMBOL': 'KEGGSYMBOL'})

    value_cnt = draw_kegg_map(pathinput,outdir,to_draw_proteomics)
    img_filename = "%s.pdf" % pathinput
                
    st.write("### KEGG Pathway with Halomonas TD1.0 genes highlighted")
    st.write("Will either show overall activity, normalized for that condition, or if subtracted case shows differences between conditions.")
    with open(os.path.join(outdir, img_filename),"rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="900" height="700" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)
    
    st.write("### Histogram for condition")
    col1, col2 = st.columns(2)
    with col1:
        PCA[input_for_conditions].hist()
        fig = plt.show()
        st.pyplot(fig)
    with col2:
        st.dataframe(value_cnt)
    
    
    
    
    
    
    
elif choice == 'Cobrapy Cell Model':
    st.markdown('# Predictions using omic-constrained TD1.0 model')
    st.markdown('Using the halomonas cobrapy model, we will constrain the upper bound of all reactions based on their protein quantity at 20, 60, and 100 g/L salt conditions.')
    st.markdown('- Note that there are 2271 total reactions in the model. 2148 are accounted for in the omics data, 32 are not annotated in the data. There are also 91 "pseudo" import/export reactions in the model which we can constrain separately from the omics data, in order to vary media conditions like glucose and measure effects.')
    
    st.sidebar.write('Exchange_Glucopyranose: glucose, Exchange_ACET: acetate, Exchange_CIT: citrate, BUTtex: butyrate,EX_SUC_e: succinate, FRUtex: fructose')
    carbon_sources = ['Exchange_Glucopyranose','Exchange_ACET','GLYCtex','Exchange_L-LACTATE','Exchange_CIT','Exchange_Trehalose','BUTtex','EX_SUC_e','ETOHtex','FRUtex']
    input_for_conditions = st.sidebar.selectbox('Select Carbon Source:', carbon_sources)
    yeast_adjust = st.sidebar.radio("Media Selection: ",('Minimal','LB','LB+'))
    
    reset_carbon(lb60_model,input_for_conditions,yeast_adjust)
    reset_carbon(lb20_model,input_for_conditions,yeast_adjust)
    reset_carbon(lb100_model,input_for_conditions,yeast_adjust)
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown('## Testing growth, PHA, and ectoine')
        s6, e6, p6 = testing_model(lb60_model)
        s2, e2, p2 = testing_model(lb20_model)
        s10, e10, p10 = testing_model(lb100_model)

        d = {'Growth': [s6.objective_value,s2.objective_value,s10.objective_value],
             'Ectoine mmol/w': [e6.objective_value,e2.objective_value,e10.objective_value],
             'PHA mmol/w': [p6.objective_value,p2.objective_value,p10.objective_value]}
        df = pd.DataFrame(data=d, index=['60 g/L', '20 g/L', '100 g/L'])
        st.write(df)
        
    with col2:
        st.markdown('## Exchange reactions at each condition')
        fluxes_df=data_frame_output_flux(lb60_model,'lb60_flux')
        fluxes_df_2=data_frame_output_flux(lb20_model,'lb20_flux')
        fluxes_df_3=data_frame_output_flux(lb100_model,'lb100_flux')
        fluxes_df=pd.concat([fluxes_df_2,fluxes_df_3,fluxes_df], axis=1)
        to_plot=fluxes_df[(fluxes_df['lb20_flux']!=0) | (fluxes_df['lb100_flux']!=0) | (fluxes_df['lb60_flux']!=0)]
        st.dataframe(to_plot.sort_values(by=['lb60_flux']), 1000, 300)
    
    col1,col2,col3 = st.columns(3)
    with col1:
        fig = testing_amino_acids(lb60_model,lb20_model,lb100_model,input_for_conditions,media=input_for_conditions,
                            lb=0,ub=20,a=2,kind='Growth',both=False)
        st.pyplot(fig.figure)
        
    with col2:
        fig_2 = testing_amino_acids(lb60_model,lb20_model,lb100_model,input_for_conditions,media=input_for_conditions,
                            lb=0,ub=20,a=2,kind='Ectoine',both=False)
        st.pyplot(fig_2.figure)
    with col3:
        fig_3 = testing_amino_acids(lb60_model,lb20_model,lb100_model,input_for_conditions,media=input_for_conditions,
                            lb=0,ub=20,a=2,kind='PHA',both=False)
        st.pyplot(fig_3.figure)
    
    st.markdown('## Predicted production of Amino Acids')
    rgs_aa = pd.DataFrame(columns = ['AA','lb60','lb20','lb100'])
    rgs_aa['AA'] = amino_acids
    rgs_aa['lb60'] = test_aa_production(lb60_model).Output
    rgs_aa['lb20'] = test_aa_production(lb20_model).Output
    rgs_aa['lb100'] = test_aa_production(lb100_model).Output
    fig_4 = rgs_aa.plot(x="AA", y=["lb20", "lb100", "lb60"], kind="bar", figsize=(19,7))
    st.pyplot(fig_4.figure)
    
    
    #st.write(lb100_model.summary())
    
    
