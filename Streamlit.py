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
from IPython.display import Image, HTML
import pandas as pd
import numpy as np
import os
import altair as alt

def to_df(result):
    return pd.read_table(io.StringIO(result),names=['KEGG Gene','Name'])
def to_df2(result):
    return pd.read_table(io.StringIO(result))

def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

# A bit of helper code to shorten long text
def head(text, lines=10):
    return text.split("\n")[:lines]

def draw_kegg_map(map_id, outdir,df):

    pathway = KGML_parser.read(kegg_get(map_id, "kgml"))
    canvas = KGMLCanvas(pathway, import_imagemap=True)
    img_filename = "%s.pdf" % map_id
    canvas.draw(os.path.join(outdir, img_filename))

#import the data
matplotlib.pyplot.switch_backend('Agg') 
read_and_cache_csv = st.cache(pd.read_csv)
countData = read_and_cache_csv('https://raw.githubusercontent.com/helloftroy/Omics_TD1.0/main/R_Landing.csv')
scatterData = read_and_cache_csv('https://raw.githubusercontent.com/helloftroy/Omics_TD1.0/main/Scatter_New.csv')
countData = pd.DataFrame(countData)
golist = list(dict.fromkeys(scatterData.GONAME))
pathlist = list(dict.fromkeys(countData.KEGGNamesofPath))


st.write("# Halomonas TD1.0 Omics")
#col1, col2, col3 = st.columns(3)
#col1.metric("GO Term", countData_filtered.iloc[1,36])
#col2.metric("Number of Genes" , countData_filtered.shape[0])
#col3.metric("Total Number Genes", countData.shape[0])
choice = st.sidebar.radio("Chart: ",('heatmap','scatterplot','MAP'))
vis_choice = st.sidebar.radio("Selection ",('NaCl 20, 60, 100','Fermentation over time'))


if vis_choice == 'NaCl 20, 60, 100':
    condition_a = 'NaCl60gLRPKM'
    condition_b = 'NaCl20gLRPKM'
    condition_c = 'NaCl100g/LRPKM'
    condition_aa = '20gNaCl'
    condition_bb = '100gNaCl'
    condition_cc = 'High_Nitrogen'
    labels = ['LB60','LB20','LB100']
    pcondition_a = '60gLProteinpgDW'
    pcondition_b = '20gLProteinpgDW'
    pcondition_c = 'NaCl100g/LProteinpgDW'
    pcondition_aa = '20gNaCl_P'
    pcondition_bb = '100gNaCl_P'
    pcondition_cc = 'High_N_P'
    vscondition_a = 'LB100vs'
    vscondition_b = 'LB20vs'
    vscondition_c = 'HighNvs'
    xlabels = 'LB60'
    ylabels = 'LB20, LB100, High N'

elif vis_choice == 'Fermentation over time':
    condition_a = 'Fermentation9hRPKM'
    condition_b = 'Fermentation19hRPKM'
    condition_c = 'Fermentation30hRPKM'
    labels = ['9hr','19hr','30hr']
    condition_aa = '19hr'
    condition_bb = '30hr'
    condition_cc = 'x'
    pcondition_a = 'Fermentation9hProteinpgDW'
    pcondition_b = 'Fermentation19hProteinpgDW'
    pcondition_c = 'Fermentation30hProteinpgDW'
    pcondition_aa = '19hr_P'
    pcondition_bb = '30hr_P'
    pcondition_cc = 'x'
    vscondition_a = '9hrvs'
    vscondition_b = '19hrvs'
    vscondition_c = '30hrvs'
    xlabels = '9 hours'
    ylabels = '19, 30 hours'


if choice == 'heatmap':
    st.write(""" ### Heatmap of different Treatments""")
    kegginput = st.selectbox('Select Pathway KEGG:', pathlist)
    #kegginput = st.text_input("Select Pathway KEGG",max_chars=200)
    truth = []
    for key in countData.KEGGPATH.fillna('a'):
        if kegginput in key and key is not float:
            truth.append(True)
        else:
            truth.append(False)
    countData_filtered = countData[truth]
    logfiltered = countData_filtered[['Name','NaCl60gLRPKM','NaCl20gLRPKM','NaCl100g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM','Fermentation30hRPKM']]
    logfiltered[['NaCl60gLRPKM','NaCl20gLRPKM','NaCl100g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM','Fermentation30hRPKM']] = np.log2(logfiltered[['NaCl60gLRPKM','NaCl20gLRPKM','NaCl100g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM','Fermentation30hRPKM']])
    ylabel = logfiltered['Name'].str[:28]
    
    exp_filtered = countData_filtered[['Name','60gLProteinpgDW','20gLProteinpgDW','NaCl100g/LProteinpgDW','Fermentation9hProteinpgDW','Fermentation19hProteinpgDW','Fermentation30hProteinpgDW']]
    exp_filtered[['60gLProteinpgDW','20gLProteinpgDW','NaCl100g/LProteinpgDW','Fermentation9hProteinpgDW','Fermentation19hProteinpgDW','Fermentation30hProteinpgDW']] = np.exp(logfiltered[['NaCl60gLRPKM','NaCl20gLRPKM','NaCl100g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM','Fermentation30hRPKM']])

    
    
    col1, col2 = st.columns(2)
    with col1:
        st.write(""" #### Transcriptomic Data""")
        sns_plot = sns.clustermap(logfiltered[[condition_a,condition_b,condition_c]],
                              xticklabels=labels, yticklabels=ylabel,method='ward',
                                  standard_scale=1,cmap="coolwarm",figsize=(5,9))
        sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize = 12)
        sns_plot.ax_heatmap.set_yticklabels(sns_plot.ax_heatmap.get_ymajorticklabels(), fontsize = 5)
        x0, _y0, _w, _h = sns_plot.cbar_pos
        sns_plot.ax_cbar.set_position([x0, 0.9, 0.04, 0.15])
        st.pyplot(sns_plot)

        st.dataframe(logfiltered[['Name',condition_a,condition_b,condition_c]].sort_values(condition_a, ascending=False).style.hide_index().highlight_max(axis=0))

    
    with col2:
        st.write(""" #### Proteomic Data""")
        
        sns_plot_p = sns.clustermap(exp_filtered[[pcondition_a,pcondition_b,pcondition_c]],
                              xticklabels=labels, yticklabels=ylabel,method='ward',
                                  standard_scale=1,cmap="coolwarm",figsize=(5,9))
        sns_plot_p.ax_heatmap.set_xticklabels(sns_plot_p.ax_heatmap.get_xmajorticklabels(), fontsize = 12)
        sns_plot_p.ax_heatmap.set_yticklabels(sns_plot_p.ax_heatmap.get_ymajorticklabels(), fontsize = 5)
        x0, _y0, _w, _h = sns_plot_p.cbar_pos
        sns_plot_p.ax_cbar.set_position([x0, 0.9, 0.04, 0.15])
        st.pyplot(sns_plot_p)
        st.dataframe(exp_filtered[['Name',pcondition_a,pcondition_b,pcondition_c]].sort_values(pcondition_a, ascending=False).style.hide_index().highlight_max(axis=0))

elif choice == 'scatterplot':
    st.write("""### Scatterplot of different Treatments""")
    col1, col2 = st.columns(2)
    with col1:
        GO = st.selectbox('Select GO TERM:', golist)
    with col2:
        st.write(""" **Hover** over points to see genes!""")  
        st.write(""" **Select** GO term at left""")      
#dna_protein = st.radio("Omics: ",('DNA','Protein','DNA vs Protein'))
    
    scatter_filtered = scatterData[scatterData.GONAME==GO]  
    X_plot=np.linspace(0,4)
    Y_plot=X_plot
    fig=plt.plot(X_plot, Y_plot, linewidth=1, color='k')
    plt.show()
    
    with col1:
  #  if dna_protein == 'DNA':
        plt.clf()
        st.write("""#### Transcriptomics""")
        scatter_filtered_d = scatter_filtered[(scatter_filtered.Treatment==condition_aa) | (scatter_filtered.Treatment==condition_bb) |(scatter_filtered.Treatment==condition_cc)]
        scattter=alt.Chart(scatter_filtered_d).mark_circle().encode(x='Control',y='Treatments',color=alt.Color('Treatment', scale=alt.Scale(scheme='category10')),tooltip='Fasta')
        polynomial_fit = [scattter.transform_regression("Control", "Treatments").mark_line()]  
        
        #could do polynomial fit here!  
        scattter=alt.layer(scattter, *polynomial_fit)           
        st.altair_chart(scattter)
    
    with col2:
        plt.clf()
        st.write("""#### Proteomics""")
        scatter_filtered_p = scatter_filtered[(scatter_filtered.Treatment==pcondition_aa) | (scatter_filtered.Treatment==pcondition_bb) | (scatter_filtered.Treatment==pcondition_cc)]
        scattter2 = alt.Chart(scatter_filtered_p).mark_circle().encode(x='Control',y='Treatments',color=alt.Color('Treatment', scale=alt.Scale(scheme='category10')),tooltip='Fasta')
        polynomial_fit = [scattter2.transform_regression("Control", "Treatments").mark_line()]  
        
        #could do polynomial fit here!  
        scattter2=alt.layer(scattter2, *polynomial_fit)           
        
        st.altair_chart(scattter2)
    st.write("""## Search KEGG Halomonas database""")
    kegginput = st.text_input("Seach KEGG Halomonas Database", value="Porin",max_chars=200,help='Type + between the words you are searching for')
    st.write(to_df(kegg_find('hel', kegginput).read()))   



elif choice == 'MAP':
    st.write("### List of Halo KEGG pathways")
    st.write(to_df(kegg_list('pathway', 'ko').read()))
    outdir = st.text_input("Please enter the pathway to you documents folder, in order to generate KEGG file",value="/Users/hellpark/Desktop/",help='Ex: /Users/hellpark/Desktop/')
    pathinput = st.text_input("To generate Halomonas map of genes in KEGG, please enter pathway of interest",value='ko00010',help='type pathway like, ko00010 after the path:')

    to_draw_proteomics = countData[["KEGGSYMBOL",pcondition_a]]
    to_draw_proteomics = to_draw_proteomics.rename(columns={pcondition_a: 'val', 'KEGGSYMBOL': 'KEGGSYMBOL'})

    str = to_df2(kegg_get(pathinput).read())
    draw_kegg_map(pathinput,outdir,to_draw_proteomics)
    img_filename = "%s.pdf" % pathinput

    st.write("### KEGG Pathway with Halomonas TD1.0 genes highlighted")
    st.write("Shows either LB60, or 9hr fermentation proteomics highlighted, depending on selection at left.")
    with open(os.path.join(outdir, img_filename),"rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="900" height="700" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)
    st.write("### List of all genes in pathway chosen!")
    st.table(str)
