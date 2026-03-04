"""
drug_likeness.py
================
Lipinski Ro5 + Veber ルールによる Drug-likeness 計算。
pipeline 全体でこの1関数を使う。
"""

from rdkit.Chem import Descriptors, rdMolDescriptors


def calculate_drug_likeness(mol) -> dict:
    """
    Lipinski Ro5 + Veber ルールを計算。

    Returns:
        {
          MW, LogP, HBD, HBA, PSA, RotBonds, Rings,
          Ro5    : bool,   # Lipinski Rule of Five
          Veber  : bool,   # Veber rules (PSA ≤ 140, RotBonds ≤ 10)
          DrugLike : bool, # Ro5 AND Veber
        }
    """
    mw   = Descriptors.ExactMolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    psa  = Descriptors.TPSA(mol)
    rot  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    ro5  = (mw <= 500) and (logp <= 5) and (hbd <= 5) and (hba <= 10)
    veber = (psa <= 140) and (rot <= 10)
    return {
        "MW": round(mw, 2), "LogP": round(logp, 2),
        "HBD": hbd, "HBA": hba, "PSA": round(psa, 1),
        "RotBonds": rot, "Rings": rings,
        "Ro5": ro5, "Veber": veber,
        "DrugLike": ro5 and veber,
    }
