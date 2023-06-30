# Contents of ~/my_app/pages/page_3.py
import streamlit as st
import scanpy as sc
import pandas as pd
from pyscenic.utils import load_motifs
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize
import io
import base64

@st.cache_resource
def upload_sc(adata):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return sc.read_h5ad(adata)

@st.cache_data
def upload_df(_df):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return load_motifs(_df)

@st.cache_data
def upload_df_rss(dataframe):
    # IMPORTANT: Cache the conversion to prevent computation on every rerun
    return pd.read_csv(dataframe,index_col=0)
##############################################################################################################

BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
COLUMN_NAME_LOGO = "MotifLogo"
COLUMN_NAME_MOTIF_ID = "MotifID"
COLUMN_NAME_TARGETS = "TargetGenes"

##############################################################################################################

def rss_plot(adata,rss):
        cats = sorted( list(set(adata.obs['cell_type'])), key=str )
        fig = plt.figure(figsize=(15, 12))
        for c,num in zip(cats, range(1,len(cats)+1)):
                x=rss.T[c]      
                ax = fig.add_subplot(3,6,num)
                plot_rss(rss, c, top_n=5, max_n=None, ax=ax)
                ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
                for t in ax.texts:
                        t.set_fontsize(12)
                ax.set_ylabel('')
                ax.set_xlabel('')
                adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )

        fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
        fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
        plt.tight_layout()
        plt.rcParams.update({
        'figure.autolayout': True,
                'figure.titlesize': 'large' ,
                'axes.labelsize': 'medium',
                'axes.titlesize':'large',
                'xtick.labelsize':'medium',
                'ytick.labelsize':'medium'
                })
        return fig

def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
    """
    :param df:
    :param base_url:
    """
    # Make sure the original dataframe is not altered.
    df = df.copy()
    
    # Add column with URLs to sequence logo.
    def create_url(motif_id):
        return "{}{}.png".format(base_url, motif_id)

    df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))
    
    # # Truncate TargetGenes.
    # def truncate(col_val):
    #     return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
    # df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))
    df = df.iloc[ [ True if x in motif_options else False for x in df.index.get_level_values('TF') ] ,:]
    df = df.reset_index()
    df = df.droplevel(level=0, axis=1)
    nombres =  ['TF','MotifID','AUC',
    'NES',
    'MotifSimilarityQvalue',
    'OrthologousIdentity',
    'Annotation',
    'Context',
    'TargetGenes',
    'RankAtMax',
    'Image']
    df.columns = nombres
    del df["TargetGenes"]

    return df
    


st.set_page_config(layout="wide")
st.set_option('deprecation.showPyplotGlobalUse', False)

##############################################################################################################

css = '''
        <style>
            .stTabs [data-baseweb="tab-list"] button [data-testid="stMarkdownContainer"] p {
            font-size:1.5rem;
            }
        </style>
        '''

st.markdown(css, unsafe_allow_html=True)  
st.markdown("# Visualize TF enrichment")


##############################################################################################################
def image_to_button(image,name):
            img = io.BytesIO()
            image.savefig(img, format='png',dpi=300,bbox_inches='tight')
            img.seek(0)
            b64 = base64.b64encode(img.read()).decode()
            htm = f'<a href="data:file/txt;base64,{b64}" download="{name}"><input type="button" value="Download Image"></a>'
            return htm


##############################################################################################################
if 'adata' not in st.session_state:
    f_adata_upload = st.file_uploader("Upload your Dataset",type="h5ad")
    if f_adata_upload is not None:
        with st.spinner("Calculating TF Enrichment ⏱️"):
                df_motifs = upload_df("reg.csv")
                rss_image_t = upload_df_rss("auc_mtx.csv")
                f_adata = upload_sc(f_adata_upload)
                with st.form("my_form"):
                        rss_image = regulon_specificity_scores( rss_image_t, f_adata.obs['cell_type'] )
                        rss_fig = rss_plot(f_adata,rss_image)
                        st.pyplot(rss_fig)
                        st.markdown(image_to_button(rss_fig,"jeje.png"),unsafe_allow_html=True)
                        st.divider()
                        tf = []
                        for a in df_motifs.index.tolist():
                                tf.append(a[0])
                        tf = set(tf)
                        motif_options = st.sidebar.multiselect("Select enriched TF",tf)
                        submitted = st.form_submit_button("Run")
                        if submitted:
                                st.data_editor(display_logos(df_motifs),
                                column_config={
                                        "Image": st.column_config.ImageColumn(
                                        "Preview Image", help="Streamlit app preview screenshots")},
                                hide_index=True)
                                st.divider()
                                reg_umap = []
                                for a in motif_options:
                                        reg_umap.append(f'Regulon({a}(+))')
                                reg_umap = sc.pl.umap(f_adata,color=reg_umap,return_fig=True)
                                st.pyplot(reg_umap)
                                st.markdown(image_to_button(reg_umap,"regulons_umap_plot.png"),unsafe_allow_html=True)

                                

                        

# else:
#         with st.spinner("Calculating TF Enrichment ⏱️"):
#                 df_motifs = upload_df("reg.csv")
#                 with st.form("my_form"):
#                         tf = []
#                         for a in df_motifs.index.tolist():
#                                 tf.append(a[0])
#                         tf = set(tf)
#                         motif_options = st.sidebar.multiselect("Select enriched TF",tf)
#                         submitted = st.form_submit_button("Run")
#                         if submitted:
#                                 st.data_editor(display_logos(df_motifs),
#                                 column_config={
#                                         "Image": st.column_config.ImageColumn(
#                                         "Preview Image", help="Streamlit app preview screenshots")},
#                                 hide_index=True)



        # # Every form must have a submit button.
        # submitted = st.form_submit_button("Submit")

        #     if submitted:
        # st.write("slider", slider_val, "checkbox", checkbox_val)

    