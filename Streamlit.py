#!/usr/bin/env python
# coding: utf-8

# In[1]:


import streamlit as st
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

def to_df(result):
    return pd.read_table(io.StringIO(result),names=['KEGG Gene','Name'])
def to_df2(result):
    return pd.read_table(io.StringIO(result))

def PDF(filename):
    return HTML('<iframe src=%s width=700 height=350></iframe>' % filename)

# A bit of helper code to shorten long text
def head(text, lines=10):
    return text.split("\n")[:lines]

def draw_kegg_map(map_id, outdir):
    """Render a local PDF of a KEGG map with the passed map ID."""
    # Get the background image first
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
choice = st.sidebar.radio("Chart: ",('heatmap','scatterplot','KEGG','MAP'))
vis_choice = st.sidebar.radio("Selection ",('NaCl 20, 60, 100','Fermentation over time'))


if vis_choice == 'NaCl 20, 60, 100':
    condition_a = 'NaCl60gLRPKM'
    condition_b = 'NaCl20gLRPKM'
    condition_c = 'NaCl100g/LRPKM'
    condition_aa = '20gNaCl'
    condition_bb = '100gNaCl'
    condition_cc = 'High_Nitrogen'
    labels = ['LB60 Transcript','LB20 Transcript','LB100 Transcript']
    
    pcondition_a = '60gLProteinpgDW'
    pcondition_b = '20gLProteinpgDW'
    pcondition_c = 'NaCl100g/LProteinpgDW'
    pcondition_aa = '20gNaCl_P'
    pcondition_bb = '100gNaCl_P'
    pcondition_cc = 'High_N_P'
    labels_p = ['LB60 Protein','LB20 Protein','LB100 Protein']

elif vis_choice == 'Fermentation over time':
    condition_a = 'Fermentation9hRPKM'
    condition_b = 'Fermentation19hRPKM'
    condition_c = 'Fermentation30hRPKM'
    labels = ['9hr t','19hr t','30hr t']
    condition_aa = '19hr'
    condition_bb = '30hr'
    condition_cc = 'x'
    pcondition_a = 'Fermentation9hProteinpgDW'
    pcondition_b = 'Fermentation19hProteinpgDW'
    pcondition_c = 'Fermentation30hProteinpgDW'
    labels_p = ['9hr p','19hr p','30hr p']
    pcondition_aa = '19hr_P'
    pcondition_bb = '30hr_P'
    pcondition_cc = 'x'


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
    logfiltered = np.log2(countData_filtered[['NaCl60gLRPKM','NaCl20gLRPKM','NaCl100g/LRPKM','Fermentation9hRPKM','Fermentation19hRPKM','Fermentation30hRPKM']])
    ylabel = countData_filtered['Name'].str[:40]
    
    st.write(""" #### Transcriptomic Data""")
    sns_plot = sns.clustermap(logfiltered[[condition_a,condition_b,condition_c]],
                          xticklabels=labels, yticklabels=ylabel,method='ward',
                              standard_scale=1,cmap="coolwarm")
    sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
    st.pyplot(sns_plot)
    
    st.write(""" #### Proteomic Data""")
    exp_filtered = np.exp(countData_filtered[[pcondition_a,pcondition_b,pcondition_c]])
    sns_plot_p = sns.clustermap(exp_filtered[[pcondition_a,pcondition_b,pcondition_c]],
                          xticklabels=labels, yticklabels=ylabel,method='ward',
                              standard_scale=1,cmap="coolwarm")
    sns_plot_p.ax_heatmap.set_xticklabels(sns_plot_p.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
    st.pyplot(sns_plot_p)
    #st.table(pathlist)


elif choice == 'scatterplot':
    st.write("""### Scatterplot of different Treatments""")
    
    GO = st.sidebar.selectbox('Select GO TERM:', golist)
    scatter_filtered = scatterData[scatterData.GONAME==GO]
    
    st.write("""#### Transcriptomics""")
    scatter_filtered_d = scatter_filtered[(scatter_filtered.Treatment==condition_aa) | (scatter_filtered.Treatment==condition_bb) | (scatter_filtered.Treatment==condition_cc)]
    scattter = sns.scatterplot(data=scatter_filtered_d,x="Control", y="NaCl20_100_3.6gLRPKM",hue="Treatment",palette= 'coolwarm')
    figurescatter = scattter.figure
    st.pyplot(figurescatter)
    
    st.write("""#### Proteomics""")
    scatter_filtered_p = scatter_filtered[(scatter_filtered.Treatment==pcondition_aa) | (scatter_filtered.Treatment==pcondition_bb) | (scatter_filtered.Treatment==pcondition_cc)]
    scattter2 = sns.scatterplot(data=scatter_filtered_p,x="Control", y="NaCl20_100_3.6gLRPKM",hue="Treatment",palette= 'coolwarm')
    figurescatter2 = scattter2.figure
    st.pyplot(figurescatter2)
    st.pyplot(plt.plot([0, 0], [3.5, 3.5], linewidth=2))
    

elif choice == 'KEGG':
    kegginput = st.text_input("Seach KEGG Halomonas Database", value="Porin",max_chars=200,help='Type + between the words you are searching for')
    st.write(to_df(kegg_find('hel', kegginput).read()))

elif choice == 'MAP':
    st.write("### List of Halo KEGG pathways")
    st.write(to_df(kegg_list('pathway', 'hel').read()))
    outdir = st.text_input("Please enter the pathway to you documents folder, in order to generate KEGG file",value="/Users/hellpark/Desktop/",help='Ex: /Users/hellpark/Desktop/')
    pathinput = st.text_input("To generate Halomonas map of genes in KEGG, please enter pathway of interest",help='type pathway like, hel00010 after the path:')
    str = to_df2(kegg_get(pathinput).read())
    draw_kegg_map(pathinput,outdir)
    img_filename = "%s.pdf" % pathinput

    st.write("### KEGG Pathway with Halomonas TD1.0 genes highlighted")
    with open(os.path.join(outdir, img_filename),"rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="900" height="700" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)
    st.write("### List of all genes in pathway chosen!")
    st.table(str)


