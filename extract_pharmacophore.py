"""
extract_pharmacophore.py
========================
重要ペプチド残基からファーマコフォア特徴点を抽出する。
各残基の側鎖から以下の特徴を取り出す:
  - HBA (水素結合受容体)
  - HBD (水素結合供与体)
  - HYD (疎水性)
  - ARO (芳香族)
  - POS (正電荷)
  - NEG (負電荷)
"""

import numpy as np
from typing import Optional
from Bio.PDB import PDBParser
import warnings
warnings.filterwarnings("ignore")


# 残基ごとの薬理活性原子定義
PHARMACOPHORE_RULES = {
    "GLY": [],
    "ALA": [("HYD", ["CB"])],
    "VAL": [("HYD", ["CB", "CG1", "CG2"])],
    "ILE": [("HYD", ["CB", "CG1", "CG2", "CD1"])],
    "LEU": [("HYD", ["CB", "CG", "CD1", "CD2"])],
    "MET": [("HYD", ["CB", "CG", "SD", "CE"])],
    "PHE": [("HYD", ["CB"]), ("ARO", ["CG", "CD1", "CD2", "CE1", "CE2", "CZ"])],
    "TRP": [("HBD", ["NE1"]), ("ARO", ["CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]),
            ("HYD", ["CB"])],
    "PRO": [("HYD", ["CB", "CG", "CD"])],
    "SER": [("HBD", ["OG"]), ("HBA", ["OG"])],
    "THR": [("HBD", ["OG1"]), ("HBA", ["OG1"]), ("HYD", ["CG2"])],
    "CYS": [("HYD", ["CB", "SG"])],
    "TYR": [("HBD", ["OH"]), ("HBA", ["OH"]), ("ARO", ["CG", "CD1", "CD2", "CE1", "CE2"]),
            ("HYD", ["CB"])],
    "ASN": [("HBD", ["ND2"]), ("HBA", ["OD1"])],
    "GLN": [("HBD", ["NE2"]), ("HBA", ["OE1"])],
    "ASP": [("NEG", ["OD1", "OD2"]), ("HBA", ["OD1", "OD2"])],
    "GLU": [("NEG", ["OE1", "OE2"]), ("HBA", ["OE1", "OE2"])],
    "LYS": [("POS", ["NZ"]), ("HBD", ["NZ"])],
    "ARG": [("POS", ["NH1", "NH2", "NE"]), ("HBD", ["NH1", "NH2"])],
    "HIS": [("ARO", ["CG", "ND1", "CD2", "CE1", "NE2"]), ("HBD", ["NE2"])],
}

FEATURE_COLORS = {
    "HBA": "red",
    "HBD": "blue",
    "HYD": "yellow",
    "ARO": "orange",
    "POS": "cyan",
    "NEG": "magenta",
}


def extract_residue_features(residue, rules: list) -> list[dict]:
    """1残基から特徴点を抽出"""
    features = []
    res_name = residue.get_resname()
    res_num  = residue.get_id()[1]
    atom_dict = {a.get_name(): a for a in residue.get_atoms()}

    for feat_type, atom_names in rules:
        coords = []
        for aname in atom_names:
            if aname in atom_dict:
                coords.append(atom_dict[aname].coord)
        if coords:
            center = np.mean(coords, axis=0)
            features.append({
                "type"    : feat_type,
                "resname" : res_name,
                "resnum"  : res_num,
                "atoms"   : atom_names,
                "coord"   : center,
                "color"   : FEATURE_COLORS.get(feat_type, "gray"),
            })
    return features


def extract_pharmacophore(model, chain_id: str,
                           key_residues: Optional[list[int]] = None) -> list[dict]:
    """
    ペプチドチェーンからファーマコフォア特徴を抽出。
    key_residues: 重要残基番号のリスト (Noneなら全残基)
    """
    features = []
    chain = model[chain_id]
    for residue in chain.get_residues():
        res_num = residue.get_id()[1]
        if key_residues and res_num not in key_residues:
            continue
        res_name = residue.get_resname()
        rules = PHARMACOPHORE_RULES.get(res_name, [])
        features.extend(extract_residue_features(residue, rules))
    return features


def print_pharmacophore(features: list[dict]):
    print("\n  ファーマコフォア特徴点:")
    print(f"  {'#':<3} {'タイプ':<6} {'残基':<8} {'座標 (x, y, z)'}")
    print("  " + "-" * 55)
    for i, f in enumerate(features):
        coord = f["coord"]
        print(f"  {i:<3} {f['type']:<6} {f['resname']}{f['resnum']:<5}"
              f"  ({coord[0]:7.2f}, {coord[1]:7.2f}, {coord[2]:7.2f})")


def save_pharmacophore_pml(features: list[dict], output_path: str):
    """PyMOL用ファーマコフォア可視化スクリプトを出力"""
    lines = ["from pymol import cmd, cgo", ""]
    for i, f in enumerate(features):
        x, y, z = f["coord"]
        color = f["color"]
        feat_label = f"{f['type']}_{f['resname']}{f['resnum']}"
        lines.append(
            f"cmd.pseudoatom('pharm_{i}', pos=[{x:.3f},{y:.3f},{z:.3f}], "
            f"label='{feat_label}')"
        )
        lines.append(f"cmd.color('{color}', 'pharm_{i}')")
        lines.append(f"cmd.show('sphere', 'pharm_{i}')")
    lines.append("cmd.set('sphere_scale', 0.5)")
    with open(output_path, "w") as fp:
        fp.write("\n".join(lines))
    print(f"  PyMOL スクリプト保存: {output_path}")


def save_pharmacophore_csv(features: list[dict], output_path: str):
    import csv
    with open(output_path, "w", newline="") as fp:
        writer = csv.DictWriter(fp, fieldnames=["index", "type", "resname", "resnum",
                                                  "x", "y", "z"])
        writer.writeheader()
        for i, f in enumerate(features):
            x, y, z = f["coord"]
            writer.writerow({"index": i, "type": f["type"],
                             "resname": f["resname"], "resnum": f["resnum"],
                             "x": round(float(x), 3), "y": round(float(y), 3),
                             "z": round(float(z), 3)})
    print(f"  ファーマコフォア CSV 保存: {output_path}")


if __name__ == "__main__":
    import sys
    from analyze_interactions import load_structure, find_contacts, score_peptide_residues

    pdb = sys.argv[1] if len(sys.argv) > 1 else "Protein_Peptide.pdb"
    model = load_structure(pdb)

    contacts = find_contacts(model, "A", "B")
    scores   = score_peptide_residues(contacts)
    # スコア上位の残基番号を取得 (スコア>0のもの)
    key_resnums = [num for (num, _), score in scores.items() if score > 0]

    features = extract_pharmacophore(model, "B", key_residues=key_resnums)
    print_pharmacophore(features)
    save_pharmacophore_csv(features, "pharmacophore.csv")
    save_pharmacophore_pml(features, "pharmacophore.pml")
