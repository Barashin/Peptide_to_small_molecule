"""
residue_defs.py
===============
アミノ酸の分類・SMILES・ファーマコフォア定義を一元管理する。
"""

# ──────────────────────────────────────────
# アミノ酸分類 (analyze_interactions / generate_pocket_molecules / pharmacophore_bridge 共通)
# ──────────────────────────────────────────

HYDROPHOBIC_AA = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "TRP", "PRO", "TYR"}
POSITIVE_AA    = {"LYS", "ARG", "HIS"}
NEGATIVE_AA    = {"ASP", "GLU"}
POLAR_AA       = {"SER", "THR", "ASN", "GLN", "CYS", "TYR", "TRP"}

# H結合供与体・受容体原子名
HBOND_DONORS    = {"N", "ND1", "ND2", "NE", "NE1", "NE2", "NH1", "NH2",
                   "NZ", "OG", "OG1", "OH"}
HBOND_ACCEPTORS = {"O", "OD1", "OD2", "OE1", "OE2", "OG", "OG1",
                   "OH", "NE2", "ND1"}

# 残基 → 機能タイプ (POS/NEG/ARO/HYD/POL)
RES_TYPE = {
    "LYS": "POS", "ARG": "POS", "HIS": "POS",
    "ASP": "NEG", "GLU": "NEG",
    "PHE": "ARO", "TYR": "ARO", "TRP": "ARO",
    "LEU": "HYD", "ILE": "HYD", "VAL": "HYD",
    "ALA": "HYD", "MET": "HYD", "PRO": "HYD",
    "SER": "POL", "THR": "POL", "ASN": "POL",
    "GLN": "POL", "CYS": "POL",
    "GLY": "HYD",
}

# ──────────────────────────────────────────
# 残基 → SMILES フラグメント (側鎖部分, 結合点は [*])
# ──────────────────────────────────────────

RESIDUE_SMILES = {
    "GLY": None,
    "ALA": "C[*]",
    "VAL": "CC(C)[*]",
    "ILE": "CCC(C)[*]",
    "LEU": "CC(C)C[*]",
    "MET": "CSCC[*]",
    "PHE": "c1ccccc1C[*]",
    "TRP": "c1ccc2[nH]ccc2c1C[*]",
    "PRO": "C1CC[NH2+]C1[*]",
    "SER": "OC[*]",
    "THR": "OC(C)[*]",
    "CYS": "SC[*]",
    "TYR": "Oc1ccc(C[*])cc1",
    "ASN": "NC(=O)C[*]",
    "GLN": "NC(=O)CC[*]",
    "ASP": "OC(=O)C[*]",
    "GLU": "OC(=O)CC[*]",
    "LYS": "NCCCC[*]",
    "ARG": "NC(=N)NCCC[*]",
    "HIS": "c1cnc[nH]1C[*]",
}

# 残基ごとのフラグメント「リーチ」推定値 (Å)
FRAG_REACH = {
    "ALA": 1.3, "VAL": 2.5, "ILE": 3.8, "LEU": 3.8,
    "MET": 3.8, "PHE": 4.5, "TRP": 5.5, "PRO": 2.5,
    "SER": 1.3, "THR": 2.5, "CYS": 2.5, "TYR": 5.0,
    "ASN": 3.8, "GLN": 5.0, "ASP": 3.8, "GLU": 5.0,
    "LYS": 5.0, "ARG": 6.3, "HIS": 4.5, "GLY": 0.0,
}
