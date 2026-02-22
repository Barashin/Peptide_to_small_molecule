"""
ligand_efficiency.py
====================
HAC (Heavy Atom Count)・LE (Ligand Efficiency)・LE グレード計算。
"""

from rdkit import Chem


def calc_hac(smiles: str):
    """SMILES から重原子数を計算。失敗時は None。"""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms() if mol else None


def calc_le(score, hac):
    """LE = |affinity| / HAC。"""
    if score is None or hac is None or hac == 0:
        return None
    return round(abs(score) / hac, 4)


def le_grade(le) -> str:
    """LE を日本語グレード文字列に変換。"""
    if le is None:
        return "N/A"
    if le >= 0.4:
        return "◎ 優秀"
    if le >= 0.3:
        return "○ 良好"
    if le >= 0.2:
        return "△ 許容"
    return "✕ 非効率"
