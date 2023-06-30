# Contents of ~/my_app/pages/page_3.py
import streamlit as st
import scanpy as sc
import matplotlib.pyplot as plt
import io
import base64


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
css = '''
        <style>
            .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {
            font-size:1.5rem;
            }
        </style>
        '''

st.markdown(css, unsafe_allow_html=True)  
st.markdown("# Differential Expression Analysis")


@st.cache_resource
def upload_sc(adata):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return sc.read_h5ad(adata)


if "ttest" not in st.session_state:
    st.session_state.ttest = "t-test"
    st.session_state.wilcoxon = False


ngenes = st.sidebar.slider("n_genes",1,10,4)
rank_method = st.sidebar.radio(options=['t-test', 'wilcoxon' ],
                    label="Rank genes for characterizing groups.",
                    help="The default method is 't-test', 'wilcoxon' uses Wilcoxon rank-sum",
                    key="ttest",
                    horizontal=True)

if 'adata' not in st.session_state:
    f_adata_upload = st.file_uploader("Upload your Dataset",type="h5ad")
    if f_adata_upload is not None:
            with st.spinner("Calculating DEA ⏱️"):
                f_adata = upload_sc(f_adata_upload)
                f_adata.raw.var.index = f_adata.var.index
                sc.pp.highly_variable_genes(f_adata,flavor="seurat",n_top_genes=2000)
                f_adata = f_adata[:, f_adata.var.highly_variable]
                f_adata.raw = f_adata
                sc.pp.normalize_total(f_adata, target_sum=1e4)
                sc.pp.log1p(f_adata)
                sc.pp.scale(f_adata, max_value=10)
                sc.tl.rank_genes_groups(f_adata,method=rank_method,groupby="cell_type")
                tab1, tab2 = st.tabs(["Mean Gene Values |", "Log Fold Changes |"])
                with tab1:
                    mgv = sc.pl.rank_genes_groups_dotplot(f_adata, n_genes=ngenes,return_fig=True)
                    st.pyplot(mgv)
                    st.markdown(image_to_button(mgv,"mean_gene_values_plot.png"),unsafe_allow_html=True)
                with tab2:
                    lfc = sc.pl.rank_genes_groups_dotplot(f_adata, n_genes=ngenes, values_to_plot='logfoldchanges', min_logfoldchange=2, vmax=3, vmin=-3, cmap='bwr',return_fig=True)
                    st.pyplot(lfc)
                    st.markdown(image_to_button(lfc,"log_fold_changes_plot"),unsafe_allow_html=True)
else:
    with st.spinner("Calculating DEA ⏱️"):
        st.session_state['adata'].raw.var.index = st.session_state['adata'].var.index
        sc.pp.highly_variable_genes(st.session_state['adata'],flavor="seurat",n_top_genes=2000)
        st.session_state['adata'] = st.session_state['adata'][:, st.session_state['adata'].var.highly_variable]
        st.session_state['adata'].raw = st.session_state['adata']
        sc.pp.normalize_total(st.session_state['adata'], target_sum=1e4)
        sc.pp.log1p(st.session_state['adata'])
        sc.pp.scale(st.session_state['adata'])
        sc.tl.rank_genes_groups(st.session_state['adata'],method=rank_method,groupby="cell_type")
        tab1, tab2 = st.tabs(["Mean Gene Values |", "Log Fold Changes |"])
        with tab1:
            mgv = sc.pl.rank_genes_groups_dotplot(st.session_state['adata'], n_genes=ngenes,return_fig=True)
            st.pyplot(mgv)
            st.markdown(image_to_button(mgv,"mean_gene_values_plot.png"),unsafe_allow_html=True)
        with tab2:
            lfc = sc.pl.rank_genes_groups_dotplot(st.session_state['adata'], n_genes=ngenes, values_to_plot='logfoldchanges', min_logfoldchange=2, vmax=3, vmin=-3, cmap='bwr',return_fig=True)
            st.pyplot(lfc)
            st.markdown(image_to_button(lfc,"log_fold_changes_plot"),unsafe_allow_html=True)
      
