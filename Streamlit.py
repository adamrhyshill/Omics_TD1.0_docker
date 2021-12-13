#!/usr/bin/env python
# coding: utf-8

# In[1]:


import warnings
import streamlit as st
import Bio
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas
from Bio.Graphics.ColorSpiral import ColorSpiral
import matplotlib
matplotlib.use('tkagg')
import matplotlib.pyplot as plt
import seaborn as sns
import pandas.util.testing as tm
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
scatterData = read_and_cache_csv('https://raw.githubusercontent.com/helloftroy/Omics_TD1.0/main/Scatter.csv')
countData = pd.DataFrame(countData)
st.write("# Halomonas TD1.0 Proteomics")
pathlist = list(dict.fromkeys(countData.KEGGNamesofPath))


#col1, col2, col3 = st.columns(3)
#col1.metric("GO Term", countData_filtered.iloc[1,36])
#col2.metric("Number of Genes" , countData_filtered.shape[0])
#col3.metric("Total Number Genes", countData.shape[0])

choice = st.sidebar.radio("Chart: ",('heatmap','scatterplot','KEGG','MAP'))
vis_choice = st.sidebar.radio("Selection ",('NaCl 20, 60, 100','Fermentation over time'))
#kegginput = st.sidebar.selectbox('Select Pathway KEGG:', pathlist)

if vis_choice == 'NaCl 20, 60, 100':
    condition_a = 'NaCl60gLRPKM'
    condition_b = 'NaCl20gLRPKM'
    condition_c = 'NaCl100g/LRPKM'
    labels = ['LB60 Transcript','LB20 Transcript','LB100 Transcript']
    pcondition_a = '60gLProteinpgDW'
    pcondition_b = '20gLProteinpgDW'
    pcondition_c = 'NaCl100g/LProteinpgDW'
    labels_p = ['LB60 Protein','LB20 Protein','LB100 Protein']

elif vis_choice == 'Fermentation over time':
    condition_a = 'Fermentation9hRPKM'
    condition_b = 'Fermentation19hRPKM'
    condition_c = 'Fermentation30hRPKM'
    labels = ['9hr t','19hr t','30hr t']
    pcondition_a = 'Fermentation9hProteinpgDW'
    pcondition_b = 'Fermentation19hProteinpgDW'
    pcondition_c = 'Fermentation30hProteinpgDW'
    labels_p = ['9hr p','19hr p','30hr p']

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
    
    sns_plot = sns.clustermap(logfiltered[[condition_a,condition_b,condition_c]],
                          xticklabels=labels, yticklabels=ylabel,method='ward',
                              standard_scale=1,cmap="coolwarm")
    sns_plot.ax_heatmap.set_xticklabels(sns_plot.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
    st.pyplot(sns_plot)
    
    exp_filtered = np.exp(countData_filtered[[pcondition_a,pcondition_b,pcondition_c]])
    sns_plot_p = sns.clustermap(exp_filtered[[pcondition_a,pcondition_b,pcondition_c]],
                          xticklabels=labels, yticklabels=ylabel,method='ward',
                              standard_scale=1,cmap="coolwarm")
    sns_plot_p.ax_heatmap.set_xticklabels(sns_plot_p.ax_heatmap.get_xmajorticklabels(), fontsize = 10)
    st.pyplot(sns_plot_p)
    #st.table(pathlist)


elif choice == 'scatterplot':
    st.write("""### Scatterplot of different Treatments""")
    golist = list(dict.fromkeys(countData.GONAME_1))
    GO = st.sidebar.selectbox('Select GO TERM:', golist)
    scatter_filtered = scatterData[scatterData.GONAME==GO]
    scattter = sns.scatterplot(data=scatter_filtered,x="NaCl60gLRPKM", y="NaCl20_100_3.6gLRPKM",hue="Treatment",palette= 'coolwarm')
    figurescatter = scattter.figure
    st.pyplot(figurescatter)
elif choice == 'KEGG':
    kegginput = st.text_input("Seach KEGG Halomonas Database", value="Porin",max_chars=200,help='Type + between the words you are searching for')
    st.write(to_df(kegg_find('hel', kegginput).read()))
    st.write(to_df(kegg_list('pathway', 'hel').read()))

elif choice == 'MAP':

    pathinput = st.text_input("Generate Halomonas Map",help='type pathway like, hel00010 after the path:')
    str = to_df2(kegg_get(pathinput).read())
    st.write(str)
    outdir = "/Users/hellpark/Desktop/"
    draw_kegg_map(pathinput,outdir)
    img_filename = "%s.pdf" % pathinput

    with open(os.path.join(outdir, img_filename),"rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="900" height="700" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)





