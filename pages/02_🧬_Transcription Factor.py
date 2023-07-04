# Contents of ~/my_app/pages/page_3.py
import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import io
import base64
import pickle 
from PIL import Image
import pandas as pd

############################################################################################################
df = pd.read_csv("todos.csv",usecols=["cell_type", "regulon","Z"])

#############################################################################################################
with open('regulons.pkl', 'rb') as f:
    regulons_list = pickle.load(f)
##############################################################################################################
def image_to_button(fig,name):
            buttons = []
            for pic in fig:
                 htm = f'<a href="data:file/txt;base64,{pic}" download="{name}"><input type="button" value="Download Image"></a>'
                 buttons.append(htm)
            return buttons

st.set_page_config(layout="wide")
st.set_option('deprecation.showPyplotGlobalUse', False)
sc.set_figure_params(dpi=300)

@st.cache_resource
def upload_sc(adata):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return sc.read_h5ad(adata)
##############################################################################################################
st.dataframe(df,use_container_width=True)
st.divider()
top = st.sidebar.multiselect("Top Regulons",regulons_list)
encoded_image_list = []
img_list = []
if top:
    for i in top:
         with io.open(f"regulon_figs/Regulon({i}(+)).png", 'rb') as archivo:
            imagen_bytes = archivo.read()
            encoded_image = base64.b64encode(imagen_bytes).decode('utf-8')
            imagen = Image.open(io.BytesIO(imagen_bytes))
            encoded_image_list.append(encoded_image)
            img_list.append(imagen) 
            # st.image(imagen,width=420)
            # st.markdown(image_to_button(encoded_image,f"regulon_{i}.png"),unsafe_allow_html=True)
st.image(img_list,width=420,caption=st.markdown(image_to_button(encoded_image_list,f"regulon_{i}.png"),unsafe_allow_html=True))
# st.markdown(image_to_button(encoded_image,f"regulon_{i}.png"),unsafe_allow_html=True)
    
   

     



