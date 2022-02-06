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
choice = st.sidebar.radio("Chart: ",('heatmap','scatterplot','KEGG','MAP', 'Model'))
vis_choice = st.sidebar.radio("Selection ",('NaCl 20, 60, 100','Fermentation over time'))



elif choice == 'MAP':
    st.write("### List of Halo KEGG pathways")

    pathinput = 'hel00010'
    img_filename = "ko00020.pdf"
    #st.write(img_filename, pathinput)
    
    pathway = KGML_parser.read(kegg_get(pathinput, "kgml"))
    canvas = KGMLCanvas(pathway, import_imagemap=True)
    
    #outdir = os.path.split('https://github.com/helloftroy/Omics_TD1.0/blob/main/')
    outdir = 'https://github.com/helloftroy/Omics_TD1.0/blob/main/'

    
    #'../omics_td1.0/outputdata/'
    #canvas.draw(os.path.join(outdir, img_filename))
    
    #draw_kegg_map(pathinput)

    st.write("### KEGG Pathway with Halomonas TD1.0 genes highlighted")
    

    folder = './outputdata'
    no_files = len(os.listdir(folder))
    file = ''
   # file_path = os.path.join(folder,str(file))
   # st.write(file_path)
   # st.write(os.path.join(file_path,os.listdir(file_path)[0]))
   # image_1 = os.path.join(file_path,os.listdir(file_path)[0])
   # st.write(image_1)
    
    with open('./outputdata/ko00020.pdf',"rb") as f:
        base64_pdf = base64.b64encode(f.read()).decode('utf-8')         
        pdf_display = F'<iframe src="data:application/pdf;base64,{base64_pdf}" width="900" height="700" type="application/pdf"></iframe>'
        st.markdown(pdf_display, unsafe_allow_html=True)
    st.write("### List of all genes in pathway chosen!")
    #st.table(str)




