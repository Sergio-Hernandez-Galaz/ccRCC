# Contents of ~/my_app/pages/page_3.py
import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import io
import base64
import pandas as pd

############################################################################################################
df = pd.read_csv("todos.csv",usecols=["cell_type", "regulon","Z"])

##############################################################################################################
renaming_dict = {0:"Trm-Exh1_VCAM1",
              1:"Tcm1_ANXA1",
              2:"Trm1_GZMH",
              3:"Trm2_TNFSF9",
              4:"Trm-Prog1_ZNF683",
              5:"Tc17_CEBPD",
              6:"Trm-Exh2_CXCL13",
              7:"Tcm2_IL7R",
              8:"Tem_GNLY",
              9:"Trm-Exh3_CCL4",
              10:"Trm-Exh4_TOX2",
              11:"Cycling-S_MCM7",
              12:"Cycling-G2/M_TOP2A",
              13:"IRG_ISG15",
              14:"Trm3_TGFBR2",
              15:"Treg_FOXP3",
              16:"Trm-Prog2_XCL1"}

##############################################################################################################
def image_to_button(image,name,format):
            img = io.BytesIO()
            image.savefig(img, format=format,dpi=300,bbox_inches='tight')
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
st.write('## scRNA-vis')
st.write("Webapp for visualization of ccRCC single cell transcripts. Created by **Sergio Hernández** 🧬")
colormap= st.sidebar.radio("UMAP Colormap",("magma","Reds","hot","rainbow","turbo","Spectral"),horizontal=True)
st.sidebar.divider()
saveas = st.sidebar.radio("Save as",("png","svg"),horizontal=True)
display_umap_celltype = st.sidebar.radio("Display Cell Type UMAP",("No","Yes"),horizontal=True)
f_adata_upload = st.file_uploader("Upload your Dataset",type="h5ad")
if f_adata_upload is not None:
    with st.form("my_form"):
        f_adata = upload_sc(f_adata_upload)
        clustering = f_adata.obs["seurat_clusters"]   
        clustering = clustering.map(renaming_dict)
        f_adata.obs["cell_type"] = clustering 

        options = st.sidebar.multiselect("Select Genes",f_adata.var_names.tolist())
        submitted = st.form_submit_button("Run")
        if submitted:
            if display_umap_celltype == "No":
                st.subheader("Selected Genes")
                sg = sc.pl.umap(f_adata,color=options,use_raw=False,cmap=colormap,return_fig=True)
                st.pyplot(sg)
                if saveas == "png":
                    st.markdown(image_to_button(sg,f"{options}.png",saveas),unsafe_allow_html=True)
                else: 
                    st.markdown(image_to_button(sg,f"{options}.svg",saveas),unsafe_allow_html=True)
                st.divider()
            else:
                cumap = sc.pl.umap(f_adata,color="cell_type",use_raw=False,cmap=colormap,return_fig=True)
                st.pyplot(cumap)
                if saveas == "png":
                    st.markdown(image_to_button(cumap,f"{options}.png",saveas),unsafe_allow_html=True)
                else: 
                    st.markdown(image_to_button(cumap,f"{options}.svg",saveas),unsafe_allow_html=True)
                st.subheader("Selected Genes")
                sg = sc.pl.umap(f_adata,color=options,use_raw=False,cmap=colormap,return_fig=True)
                st.pyplot(sg)
                if saveas == "png":
                    st.markdown(image_to_button(sg,f"{options}.png",saveas),unsafe_allow_html=True)
                else: 
                    st.markdown(image_to_button(sg,f"{options}.svg",saveas),unsafe_allow_html=True)
                st.divider()
                 

     



