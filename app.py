from datetime import date
import pandas as pd

from chebi_tools.ChEBIStandardizer import ChEBIStandardizer
#from chebi_tools.ChEBIGraph import ChEBIGraph

import streamlit as st

st.set_page_config(layout="wide")


def today():
    return date.today().strftime("%y%m%d")

STD = ChEBIStandardizer()

st.title('Metabolite-Standardizer')

text = st.text_area('Paste metabolite names (one name per line)', value="Mevalonate\nAspartate", height=None, max_chars=None, 
   key='test', help=None, on_change=None, args=None, kwargs=None, 
   placeholder=None)
   
@st.cache_data
def convert_df(df):
   return df.to_csv(index=False).encode('utf-8')   

if st.button('Convert', key=None, help=None, on_click=None, args=None, 
             kwargs=None, type="secondary", disabled=False):
    
    tokens = text.split('\n')

    tokens = [e.strip() if e != '']
    
    df = pd.DataFrame()
    
    df = STD.process_many(tokens)

    st.dataframe(df)

    csv = convert_df(df)

    st.download_button(
        "Download", 
        csv, 
        f"{today()}-standardized-metabolites.csv", 
        "text/csv", key='download-csv'
    )