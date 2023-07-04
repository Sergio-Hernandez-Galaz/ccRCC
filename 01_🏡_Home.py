# Contents of ~/my_app/pages/page_3.py
import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import io
import base64

##############################################################################################################
def image_to_button(image,name):
            img = io.BytesIO()
            image.savefig(img, format='png',dpi=300,bbox_inches='tight')
            img.seek(0)
            b64 = base64.b64encode(img.read()).decode()
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
colormap= st.sidebar.radio("UMAP Colormap",("magma","Reds","hot","rainbow","turbo"))
f_adata_upload = st.file_uploader("Upload your Dataset",type="h5ad")
if f_adata_upload is not None:
    with st.form("my_form"):
        f_adata = upload_sc(f_adata_upload)
        options = st.sidebar.multiselect("Select Genes",f_adata.var_names.tolist())
        submitted = st.form_submit_button("Run")
        if submitted:
            st.subheader("Selected Genes")
            sg = sc.pl.umap(f_adata,color=options,use_raw=False,cmap=colormap,return_fig=True)
            st.pyplot(sg)
            st.markdown(image_to_button(sg,f"{options}.png"),unsafe_allow_html=True)
            st.divider()

     



