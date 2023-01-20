import logging
import Levenshtein

import pandas as pd
import numpy as np

from pathlib import Path as P
from rdkit import Chem

from .ChEBIDownloader import ChEBIDownloader
from .tools import standardize_smiles, get_mol_props



class ChEBI():
    def __init__(self, download_dir=None):
        self.names = None
        self._names = None
        self.smiles = None
        self.inchi_keys = None
        self.downloader = ChEBIDownloader(download_dir=download_dir)
        self.downloader.download_missing()
        self.download_dir = self.downloader.download_dir
        self.files = self.downloader.files
        
        self.fn_smiles = f'{self.download_dir}/smiles.parquet'
        self.fn_inchi_keys = f'{self.download_dir}/inchi_keys.parquet'
        
        # Prepare data
        self.prepare_data()
        
        # Load data
        self.load_data()
        
        
    def prepare_data(self):
        self.prepare_smiles()
        self.prepare_inchi_keys()
        
        
    def load_data(self):
        self.load_names()
        self.load_smiles()
        self.load_inchi_keys()
        
        
    def load_names(self): 
        #names = pd.read_csv(f'{self.download_dir}/names.tsv', sep='\t', na_filter=False)
        
        names_from_compounds = pd.read_parquet(self.files['compounds'])[['CHEBI_ACCESSION', 'SOURCE', 'NAME']].replace('null', None).dropna()
        names_from_compounds['COMPOUND_ID'] = names_from_compounds.CHEBI_ACCESSION.str.replace('CHEBI:', '').astype(int)
        names_from_compounds = names_from_compounds.drop('CHEBI_ACCESSION', axis=1)
        
        names_from_names = pd.read_parquet(self.files['names'])[['SOURCE', 'COMPOUND_ID', 'NAME', 'TYPE']]
        
        names = pd.concat([names_from_names, names_from_compounds]).drop_duplicates()
        names.index = names.NAME.str.lower()
        
        r = names[names.NAME.fillna('').str.startswith('(R)-')].copy()
        r['NAME'] = r['NAME'].fillna('').str.replace('(R)-', '', regex=False)
        names = pd.concat([names, r])
        
        names = names[~names.SOURCE.isin(['Chemical Ontology'])]
        
        names.index = names['NAME'].str.lower()
        self.names = names
        self._names = self.names[~self.names.SOURCE.isin(['CBN', 'PDBeChem'])].set_index('COMPOUND_ID')
    
    
    def load_smiles(self):
        self.smiles = pd.read_parquet(self.fn_smiles)
    
    
    def load_inchi_keys(self):
        self.inchi_keys = pd.read_parquet(self.fn_inchi_keys)
    
    
    def prepare_smiles(self):
        if P(self.fn_smiles).is_file():
            return
        logging.warning('Preparing smiles')
        structures = pd.read_parquet(self.files['structures'])
        smiles = structures[structures.TYPE == 'SMILES'].set_index('COMPOUND_ID')[['STRUCTURE']]
        smiles.to_parquet(self.fn_smiles)
        
        
    def prepare_inchi_keys(self):
        if P(self.fn_inchi_keys).is_file():
            return
        logging.warning('Preparing inchi_keys')
        structures = pd.read_parquet(self.files['structures'])
        smiles = structures[structures.TYPE == 'InChIKey'][['STRUCTURE', 'COMPOUND_ID']].reset_index(drop=True)
        smiles.to_parquet(self.fn_inchi_keys)
           
            
    def get(self, token):
        token = self._process_token(token)

        match token:
            case int():
                return self._names.loc[token]
            case str():
                return self.names.loc[token.lower()]
                   
    def get_exact_compound_id(self, compound_id:int):
        return 
    
    def contains(self, name:str):
        return self.names[self.names.index.str.contains(name.lower())]
        
        
    def search(self, name:str):
        candidates = self.names.index
        name = name.lower()
        distances = [self.distance(candidate, name) for candidate in candidates]
        ndx = np.argmin(distances)
        return self.names.iloc[ndx]

    
    def distance(self, string_a, string_b):
        return Levenshtein.distance(string_a, string_b)
    
    
    def get_compound_id(self, name:str):
        results = None
        try:
            results = self.get(name)
            exact_match = True
        except KeyError as e:
            results = self.search(name)
            exact_match = False
            logging.warning(f'No excact match for {name}.')

        if isinstance(results, pd.core.series.Series):
            name = results['NAME']
            compound_id = results['COMPOUND_ID']
            return int(compound_id), exact_match
        else:
            return int(results['COMPOUND_ID'].value_counts().index[0]), exact_match
    
    
    def get_names_of_compound_id(self, compound_id:int):
        return self._names.loc[compound_id]
    
    def get_compound_id_of_name(self, name:str):
        return self.names.loc[name, 'COMPOUND_ID']
    
    def get_smiles(self, token):
        if isinstance(token, str):
             compound_id, exact_match = self.get_compound_id(token)
        else:
            compound_id = token
        try:
            return self.smiles.loc[compound_id].STRUCTURE
        except KeyError:
            return None
    
    
    def inchi_key_to_compound_id(self, inchi_key):
        return self.inchi_keys.set_index('STRUCTURE').loc[inchi_key]
    
    
    def get_synonyms(self, compound_id:int) -> list:
        try:
            synonyms = self._names.loc[compound_id].NAME
        except KeyError:
            synonyms = self.names[self.names.COMPOUND_ID == compound_id].NAME
        if isinstance(synonyms, str):
            synonyms = [synonyms]
        else:
            synonyms = list(synonyms)
        return synonyms
    
    
    def best_name(self, compound_id:int, best_length=10):
        synonyms = self.get_all_names(compound_id)
        # if just one name, return it
        if len(synonyms) == 1:
            return synonyms[0]
        ndx_alpha = [ e.replace(' ', '').isalpha() for e in synonyms ]
        if sum(ndx_alpha) > 0:
            synonyms = [ e for e in synonyms if e.replace(' ', '').isalpha()]
        ndx_shortest = np.argmin([np.abs(len(e)-best_length) for e in synonyms])
        shortest_name = synonyms[ndx_shortest]
        return shortest_name
    
    def standardize_name(self, name:str):
        compound_id, exact_match = self.get_compound_id(name) 
        new_name =  self.best_name(compound_id)
        if new_name.replace(' ','').isalpha() and not new_name.isupper():
            new_name = new_name.capitalize()
        return new_name
    
    def get_all_names(self, compound_id:int):
        return list(set(self._names.loc[compound_id, 'NAME']))
        
    def get_data(self, token):
        token = self._process_token(token)

        match token:
            case int():
                compound_id = int(token)
                exact_match = True
                name = self.best_name(compound_id)
            case str():
                name = token 
                compound_id, exact_match = self.get_compound_id(name)
        
        standardized_name = self.standardize_name(name)        
        smiles = self.get_smiles(compound_id)
        smiles_std = standardize_smiles(smiles)
        results = {
            'query_name': name,
            'standardized_name': standardized_name,
            'smiles': smiles,
            'standardized_smiles': smiles_std,
            'ChEBI': compound_id,
            'exact_match': exact_match
        }
        if smiles is not None:
            mol = self.get_mol_from_smiles(smiles_std)
            properties = self.get_mol_props(mol)
            results.update(properties)
        return results

    
    def get_molecule(self, token):
        token = self._process_token(token)
            
        match token:
            case int():
                smiles = self.get_smiles(token)
            case str():
                compound_id, exact_match = self.get_compound_id(token)
                smiles = self.get_smiles(compound_id)
                
        return self.get_mol_from_smiles(smiles)
    
    
    def get_mol_from_smiles(self, smiles: str):
        return Chem.MolFromSmiles(smiles)
    
    def get_mol_props(self, mol):
        return get_mol_props(mol)

    def get_all_compound_ids(self, compound_id):
        synonyms = self.get_all_names(compound_id)
        return self._names.loc[self._names.NAME.isin(synonyms)].index.unique().to_list()
    
    def get_best_compound_id(self, compound_id):
        return self.get_all_compound_ids(compound_id)[0]
    
    def _process_token(self, token):
        if isinstance(token, str):
            token = token.lower().replace('chebi:', '')
            if token.isnumeric():
                token = int(token)
        return token
    
    def token_to_name_and_compound_id(self, token):
        token = self._process_token(token)

        name, compound_id = None, None        
        match token:
            case int():
                compound_id = self.get_best_compound_id(token)
                name = self.best_name(compound_id)
            case str():
                name = token
                compound_id, exact_match = self.get_compound_id(name)

        return name, compound_id 
