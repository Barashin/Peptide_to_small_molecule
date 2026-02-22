"""
analyze_interactions.py
=======================
タンパク質-ペプチド複合体の相互作用を解析するモジュール。
- H結合 (N/O-H...N/O, 距離 ≤3.5Å, 角度 ≥120°)
- 疎水性接触 (C...C, 距離 ≤4.5Å)
- 静電相互作用 (charged residues, 距離 ≤6.0Å)
"""

import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
from Bio.PDB.vectors import calc_angle
import warnings
warnings.filterwarnings("ignore")


from utils.residue_defs import (HYDROPHOBIC_AA, POSITIVE_AA, NEGATIVE_AA,
                                POLAR_AA, HBOND_DONORS, HBOND_ACCEPTORS)


def load_structure(pdb_path: str, model_id: int = 0):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    return structure[model_id]


def get_chain_atoms(model, chain_id: str) -> list:
    """指定チェーンの全重原子を取得"""
    atoms = []
    for residue in model[chain_id].get_residues():
        for atom in residue.get_atoms():
            if atom.element != "H":
                atoms.append(atom)
    return atoms


def find_contacts(model, chain_protein: str, chain_peptide: str,
                  cutoff: float = 4.5) -> list[dict]:
    """
    タンパク質-ペプチド間の重原子接触を検索。
    Returns: list of contact dicts
    """
    all_atoms = list(model.get_atoms())
    ns = NeighborSearch(all_atoms)

    peptide_atoms = get_chain_atoms(model, chain_peptide)
    contacts = []

    for p_atom in peptide_atoms:
        p_res = p_atom.get_parent()
        nearby = ns.search(p_atom.coord, cutoff, level="A")
        for n_atom in nearby:
            n_res = n_atom.get_parent()
            # 異なるチェーン間のみ
            if n_res.get_parent().id != chain_protein:
                continue
            if n_atom.element == "H":
                continue
            dist = np.linalg.norm(p_atom.coord - n_atom.coord)
            contacts.append({
                "peptide_res"    : p_res.get_resname(),
                "peptide_resnum" : p_res.get_id()[1],
                "peptide_atom"   : p_atom.get_name(),
                "protein_res"    : n_res.get_resname(),
                "protein_resnum" : n_res.get_id()[1],
                "protein_atom"   : n_atom.get_name(),
                "distance"       : round(dist, 3),
                "peptide_coord"  : p_atom.coord.copy(),
                "protein_coord"  : n_atom.coord.copy(),
            })

    return contacts


def classify_contact(contact: dict) -> str:
    """接触をタイプ別に分類"""
    pa = contact["peptide_atom"]
    na = contact["protein_atom"]
    pr = contact["peptide_res"]
    d  = contact["distance"]

    if pa[0] in ("N", "O") and na[0] in ("N", "O") and d <= 3.5:
        return "hbond"
    if pa[0] == "C" and na[0] == "C" and d <= 4.5:
        if pr in HYDROPHOBIC_AA:
            return "hydrophobic"
        return "van_der_waals"
    if pr in POSITIVE_AA or pr in NEGATIVE_AA:
        if contact["protein_res"] in POSITIVE_AA or contact["protein_res"] in NEGATIVE_AA:
            if d <= 6.0:
                return "electrostatic"
    return "other"


def score_peptide_residues(contacts: list[dict]) -> dict:
    """
    各ペプチド残基の重要度スコアを計算。
    H結合: 3点, 疎水性: 1点, 静電: 2点, その他: 0.5点
    """
    weights = {"hbond": 3.0, "electrostatic": 2.0, "hydrophobic": 1.0,
               "van_der_waals": 0.5, "other": 0.3}
    scores = {}
    for c in contacts:
        key = (c["peptide_resnum"], c["peptide_res"])
        ctype = classify_contact(c)
        scores[key] = scores.get(key, 0.0) + weights.get(ctype, 0.0)
    return dict(sorted(scores.items(), key=lambda x: -x[1]))


def summarize(contacts: list[dict]) -> dict:
    """接触サマリーを生成"""
    typed = [classify_contact(c) for c in contacts]
    from collections import Counter
    type_counts = Counter(typed)
    residue_scores = score_peptide_residues(contacts)

    print("=" * 60)
    print("  タンパク質-ペプチド相互作用 サマリー")
    print("=" * 60)
    print(f"  総接触数     : {len(contacts)}")
    for t, n in type_counts.items():
        print(f"  {t:<20}: {n}")
    print()
    print("  ペプチド残基スコア (高いほど重要):")
    for (num, name), score in list(residue_scores.items())[:10]:
        print(f"    残基 {num:>3} {name} : {score:.1f}")
    print("=" * 60)

    return {"contacts": contacts, "type_counts": type_counts,
            "residue_scores": residue_scores}


if __name__ == "__main__":
    import sys
    pdb = sys.argv[1] if len(sys.argv) > 1 else "Protein_Peptide.pdb"
    model = load_structure(pdb)
    contacts = find_contacts(model, chain_protein="A", chain_peptide="B")
    result = summarize(contacts)
