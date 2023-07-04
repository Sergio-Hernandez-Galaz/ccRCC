# Contents of ~/my_app/pages/page_3.py
import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import io
import base64
import pickle 
from PIL import Image

#############################################################################################################
with open('regulons.pkl', 'rb') as f:
    regulons_list = pickle.load(f)
##############################################################################################################
def image_to_button(fig,name):
            b64 = base64.b64encode(fig.getvalue()).decode()
            htm = f'<a href="data:file/txt;base64,{b64}" download="{name}"><input type="button" value="Download Image"></a>'
            return htm

st.set_page_config(layout="wide")
st.set_option('deprecation.showPyplotGlobalUse', False)
sc.set_figure_params(dpi=300)

@st.cache_resource
def upload_sc(adata):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return sc.read_h5ad(adata)
##############################################################################################################
st.sidebar.write('## scRNA-vis Beta')
st.sidebar.write("Webapp for visualization of ccRCC single cell transcripts. Created by **Sergio Hernández** 🧬")

top = st.sidebar.multiselect("Top Regulons",regulons_list)
img_list = []
if top:
    for i in top:
         with io.open(f"regulon_figs/Regulon({i}(+)).png", 'rb') as archivo:
            imagen_bytes = archivo.read()
            imagen = Image.open(io.BytesIO(imagen_bytes))
            st.markdown(image_to_button(imagen,"regulon_{i}.png"),unsafe_allow_html=True)
            img_list.append(imagen)
            
    st.image(img_list,width=435)
   

     



