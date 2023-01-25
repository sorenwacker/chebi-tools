import Levenshtein

import numpy as np
import pandas as pd

from .ChEBIDownloader import ChEBIDownloader   

class ChEBIStandardizer():
    
    def __init__(self, download_dir=None):
        self.names = None
        self.reference_chebi = None
        
        self.downloader = ChEBIDownloader(download_dir=download_dir)
        self.downloader.download_missing()
        self.download_dir = self.downloader.download_dir
        self.files = self.downloader.files
                        
        # Load data
        self.load_data()    
        
        
    def load_data(self):
        self.load_names()
        self.load_reference_data()

    
    def load_names(self):         
        # To get all synonyms, the names from names.tsv and structures.tsv need to be combined
        names_from_compounds = pd.read_parquet(self.files['compounds'])[['CHEBI_ACCESSION', 'SOURCE', 'NAME']].replace('null', None).dropna()
        names_from_compounds['COMPOUND_ID'] = names_from_compounds.CHEBI_ACCESSION.str.replace('CHEBI:', '').astype(int)
        names_from_compounds = names_from_compounds.drop('CHEBI_ACCESSION', axis=1)
        names_from_names = pd.read_parquet(self.files['names'])[['SOURCE', 'COMPOUND_ID', 'NAME', 'TYPE']]
        names = pd.concat([names_from_names, names_from_compounds]).drop_duplicates()
        
        # Experimental, ma
        #r = names[names.NAME.fillna('').str.startswith('(R)-')].copy()
        #r['NAME'] = r['NAME'].fillna('').str.replace('(R)-', '', regex=False)
        #names = pd.concat([names, r])
        
        # Remove non-compound names
        names = names[~names.SOURCE.isin(['Chemical Ontology'])]
        
        names['NAME_lowercase'] = names['NAME'].str.lower()
        self.names = names
        
        
    def load_reference_data(self):
        reference_chebi = pd.read_parquet('/home/swacker/Projects/230119-sw__chebi-tools/_reference-chebis.parquet').dropna().drop_duplicates()
        reference_chebi.columns = ['ChEBI', 'NAME']
        reference_chebi['ChEBI'] = reference_chebi['ChEBI'].str.replace('CHEBI', 'ChEBI')
        reference_chebi.index = reference_chebi['ChEBI'].apply(lambda x: int(x.lower().replace('chebi:', '')))
        self.reference_chebi = reference_chebi
        
        
    def get(self, token):
        token = self._process_token(token)

        match token:
            case int() | np.int8() | np.int32() | np.int64() | np.uint8() | np.uint32() | np.uint64():
                return self._names.loc[token]
            case str():
                return self.names.loc[token.lower()]        
        
        
    def get(self, token):
        token = self._process_token(token)

        match token:
            case int() | np.int8() | np.int32() | np.int64() | np.uint8() | np.uint32() | np.uint64():
                return self._names.loc[token]
            case str():
                return self.names.loc[token.lower()]
            
            
    def get_by_id(self, compound_id:int):
        return self.names[self.names.COMPOUND_ID == compound_id]
        
        
    def get_by_name(self, name:str):
        return self.names[self.names.NAME.str.lower() == name.lower()]
    
    
    def search(self, name:str):
        if name.lower() in self.names.NAME_lowercase:
            exact_match = True
            return self.names.set_index('c').loc[name.lower()], exact_match
        else:
            exact_match = False
            candidates = self.names.NAME_lowercase
            name = name.lower()
            distances = [self.distance(candidate, name) for candidate in candidates]
            ndx = np.argmin(distances)
            return self.names.iloc[ndx], exact_match

    
    def distance(self, string_a, string_b):
        return Levenshtein.distance(string_a, string_b)
    
    
    def get_chebi_reference(self, token):
        token = self._process_token(token)
        
        match token:
            case int() | np.int8() | np.int32() | np.int64() | np.uint8() | np.uint32() | np.uint64():
                compound_id = token
                exact_match = True
            case str():
                search_result, exact_match = self.search(token)
                compound_id = search_result.COMPOUND_ID

        result = self.reference_chebi.loc[compound_id].to_dict()
        result.update(dict(exact_match=exact_match))
        return result

    def suggest_name(self, token):
        record = self.get_chebi_reference(token)
        name = record['NAME']
        chebi = record['ChEBI']
        if len(name) > 30:
            return chebi
        return name
        
    def _process_token(self, token):
        if isinstance(token, str):
            token = token.lower().replace('chebi:', '')
            if token.isnumeric():
                token = int(token)
        return token    