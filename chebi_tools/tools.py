import logging
import traceback

from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize

from rdkit import RDLogger
from rdkit.Chem import Descriptors
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# from rdkit.Chem.rdPartialCharges import ComputeGasteigerCharges


RDLogger.DisableLog("rdApp.*")


def get_mol_props(mol):
    properties = dict(
        molecular_weight=Descriptors.ExactMolWt(mol),
        formula=CalcMolFormula(mol),
        formal_charge=Chem.GetFormalCharge(mol),
    )
    return properties


def standardize_smiles(smiles):
    if (smiles == "") or (smiles is None):
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        clean_mol = rdMolStandardize.Cleanup(mol)
        parent_clean_mol = rdMolStandardize.FragmentParent(clean_mol)
        uncharger = (
            rdMolStandardize.Uncharger()
        )  # annoying, but necessary as no convenience method exists
        uncharged_parent_clean_mol = uncharger.uncharge(parent_clean_mol)
        te = rdMolStandardize.TautomerEnumerator()  # idem
        taut_uncharged_parent_clean_mol = te.Canonicalize(uncharged_parent_clean_mol)
        return Chem.MolToSmiles(taut_uncharged_parent_clean_mol)
    except Exception:
        logging.error(traceback.format_exc())
        return None
