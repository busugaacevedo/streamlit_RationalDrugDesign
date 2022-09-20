"""
Docker Container: https://hub.docker.com/r/continuumio/anaconda3
RDKit Installation: https://www.rdkit.org/docs/Install.html
"""
#import mols2grid
import pandas as pd
import streamlit as st
import streamlit.components.v1 as components
from rdkit import Chem
from rdkit.Chem.Descriptors import ExactMolWt, MolLogP, NumHDonors, NumHAcceptors

st.title("Rational design of Peptides and small ligands")

st.markdown("""
- App modified for peptides by [Brandon Usuga](https://github.com/busugaacevedo/peptides_BRANDON)
- App modified by [Chanin Nantasenamat](http://medium.dataprofessor.org) (aka [Data Professor](http://youtube.com/dataprofessor))
- Original app by [Justin Chavez](https://blog.reverielabs.com/building-web-applications-from-python-scripts-with-streamlit/)
""")

#@st.cache(allow_output_mutation=True)
#def download_dataset():
#    """Loads once then cached for subsequent runs"""

# Calculate descriptors
def calc_mw(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the molecular weight"""
    mol = Chem.MolFromSmiles(smiles_string)
    return ExactMolWt(mol)

def calc_logp(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the LogP"""
    mol = Chem.MolFromSmiles(smiles_string)
    return MolLogP(mol)

def calc_NumHDonors(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHDonors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHDonors(mol)

def calc_NumHAcceptors(smiles_string):
    """Given a smiles string (ex. C1CCCCC1), calculate and return the NumHAcceptors"""
    mol = Chem.MolFromSmiles(smiles_string)
    return NumHAcceptors(mol)




# Sidebar panel     -- Define variables --
st.sidebar.header('Set parameters')
st.sidebar.write('*Note: Display compounds having values less than the following thresholds*')

length_peptide = st.sidebar.slider(
        label="Peptide length",
        min_value=1,
        max_value=4,
        value=2,
        step=1,
)
weight_cutoff = st.sidebar.slider(
    label="Molecular weight",
    min_value=0,
    max_value=1000,
    value=500,
    step=10,
)
logp_cutoff = st.sidebar.slider(
    label="LogP",
    min_value=-10,
    max_value=10,
    value=5,
    step=1,
)
NumHDonors_cutoff = st.sidebar.slider(
    label="NumHDonors",
    min_value=0,
    max_value=15,
    value=5,
    step=1,
)
NumHAcceptors_cutoff = st.sidebar.slider(
    label="NumHAcceptors",
    min_value=0,
    max_value=20,
    value=10,
    step=1,
)
## -- END of the variables -- ##

## -- Load the database -- ##
#length_peptide=1
if length_peptide == 1:
    df = pd.read_csv("https://raw.githubusercontent.com/busugaacevedo/streamlit_RationalDrugDesign/main/1_peps.txt", sep="\t").dropna()
elif length_peptide == 2:
    df = pd.read_csv("https://raw.githubusercontent.com/busugaacevedo/streamlit_RationalDrugDesign/main/2_peps.txt", sep="\t").dropna()
elif length_peptide == 3:
    df = pd.read_csv("https://raw.githubusercontent.com/busugaacevedo/streamlit_RationalDrugDesign/main/3_peps.txt", sep="\t").dropna()
else:
    df = pd.read_csv("https://raw.githubusercontent.com/busugaacevedo/streamlit_RationalDrugDesign/main/4_peps.txt", sep="\t").dropna()
            
        #"TABLA1.txt", sep="\t"
        #"https://raw.githubusercontent.com/busugaacevedo/stream_drugs/main/f_2.txt?token=GHSAT0AAAAAABY6XZ2LRFFXZKIA6JVOKSTEYZIX5FA", sep=","
        #"https://www.cureffi.org/wp-content/uploads/2013/10/drugs.txt", sep="\t"
#    return df

## -- END of the database -- ##

## -- Computes -- ##
# Copy the dataset so any changes are not applied to the original cached version
#df = download_dataset().copy()
df["MW"] = df.apply(lambda x: calc_mw(x["smiles"]), axis=1)
df["LogP"] = df.apply(lambda x: calc_logp(x["smiles"]), axis=1)
df["NumHDonors"] = df.apply(lambda x: calc_NumHDonors(x["smiles"]), axis=1)
df["NumHAcceptors"] = df.apply(lambda x: calc_NumHAcceptors(x["smiles"]), axis=1)

## -- END of the compute -- ##

df_result = df[df["MW"] < weight_cutoff]
df_result2 = df_result[df_result["LogP"] < logp_cutoff]
df_result3 = df_result2[df_result2["NumHDonors"] < NumHDonors_cutoff]
df_result4 = df_result3[df_result3["NumHAcceptors"] < NumHAcceptors_cutoff]

st.write(df_result4.shape)
st.write(df_result4)


#raw_html = mols2grid.display(df_result4,
#                            #subset=["Name", "img"],
#                            subset=["img", "Sequence", "MW", "LogP", "NumHDonors", "NumHAcceptors"],
#                            mapping={"smiles": "SMILES", "generic_name": "Name"})._repr_html_()
#components.html(raw_html, width=900, height=1100, scrolling=False)
