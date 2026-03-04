"""
generate_pocket_molecules.py
============================
タンパク質-ペプチド界面の「タンパク質側」接触残基を3手法で統合スコアリングし、
ポケットの物理化学的性質に相補的な低分子リガンドを設計するモジュール。

設計の論理:
  タンパク質残基タイプ → 相補的リガンド特徴
  ─────────────────────────────────────────
  LYS/ARG (正電荷, HBD)  → 負電荷基 (COOH, SO₃H, PO₃H₂, テトラゾール)
  ASP/GLU (負電荷, HBA)  → 正電荷基 (アミン, グアニジン)
  LEU/ILE/VAL (疎水性)   → 疎水性コア (芳香環, シクロヘキシル, アルキル)
  TYR/PHE/TRP (芳香族)   → π-π スタッキング (芳香環)
  SER/THR/ASN (極性)     → H結合基 (OH, C=O, NH₂)
  PRO (疎水性環)          → 疎水性基

フラグメント組み立て戦略:
  Anchor  : ポケットの主要な残基タイプと結合するコア基
  Scaffold: 疎水性コア/リンカー
  Decorator: 追加の相互作用基

出力: candidate_pocket_ligands.sdf (ドッキング用)
"""

import os
import json
import csv
import numpy as np
from collections import defaultdict, Counter
from typing import Optional
from Bio.PDB import PDBParser
import warnings
warnings.filterwarnings("ignore")

from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import BRICS


# ──────────────────────────────────────────
# 残基タイプ定義
# ──────────────────────────────────────────
from utils.residue_defs import RES_TYPE

# タンパク質残基タイプ → 相補的リガンドフラグメント SMILES
COMPLEMENTARY_FRAGMENTS = {
    "POS": {   # 正電荷残基 (LYS/ARG) → アニオン性基
        "carboxylic_acid" : "OC(=O)",
        "sulfonamide"     : "S(=O)(=O)N",
        "tetrazole"       : "c1nnn[nH]1",
        "phosphonate"     : "P(=O)(O)O",
        "hydroxamic_acid" : "C(=O)NO",
    },
    "NEG": {   # 負電荷残基 (ASP/GLU) → カチオン性基
        "primary_amine"   : "N",
        "guanidinium"     : "NC(=N)N",
        "amidine"         : "NC(=N)",
    },
    "ARO": {   # 芳香族残基 (TYR/PHE/TRP) → π-π スタッキング
        "phenyl"          : "c1ccccc1",
        "pyridine"        : "c1ccncc1",
        "indole"          : "c1ccc2[nH]ccc2c1",
        "naphthalene"     : "c1ccc2ccccc2c1",
        "fluorophenyl"    : "c1ccc(F)cc1",
    },
    "HYD": {   # 疎水性残基 (LEU/ILE/VAL) → 疎水性基
        "cyclohexyl"      : "C1CCCCC1",
        "isobutyl"        : "CC(C)C",
        "isopropyl"       : "CC(C)",
        "tert_butyl"      : "C(C)(C)C",
        "methyl_phenyl"   : "Cc1ccccc1",
    },
    "POL": {   # 極性残基 (SER/THR) → H結合基
        "hydroxyl"        : "O",
        "amide"           : "C(=O)N",
        "hydroxymethyl"   : "CO",
    },
}

# テンプレート分子 (ポケット特性に合わせた既知フラグメント)
# フォーマット: {"name": str, "smiles": str, "rationale": str}
MOLECULE_TEMPLATES = [
    # ─── POS-rich pocket (LYS/ARG) + HYD + ARO ───
    # 1. 芳香族カルボン酸系 (ベンゼン環でARO/HYD残基、COOHでLYS/ARGに結合)
    {"name": "benzoic_acid",
     "smiles": "OC(=O)c1ccccc1",
     "rationale": "Phenyl→LEU/ILE(HYD), COOH→LYS/ARG(POS)"},
    {"name": "phenylacetic_acid",
     "smiles": "OC(=O)Cc1ccccc1",
     "rationale": "Phenyl→LEU/ILE/TYR(ARO), CH2-COOH→LYS/ARG(POS)"},
    {"name": "3_phenylpropionic_acid",
     "smiles": "OC(=O)CCc1ccccc1",
     "rationale": "Phenyl→HYD pocket, chain-COOH→POS残基"},
    {"name": "4_phenylbutyric_acid",
     "smiles": "OC(=O)CCCc1ccccc1",
     "rationale": "Extended hydrophobic chain+COOH, flexible linker"},

    # 2. インドール系 (TRPに類似したコア + アニオン性基)
    {"name": "indole_3_acetic_acid",
     "smiles": "OC(=O)Cc1c[nH]c2ccccc12",
     "rationale": "Indole→TYR/PHE(ARO)+HYD, IAA backbone→LYS/ARG(POS)"},
    {"name": "indole_3_propionic_acid",
     "smiles": "OC(=O)CCc1c[nH]c2ccccc12",
     "rationale": "Indole core+propionate, good size for pocket"},
    {"name": "5_fluoroindole_3_acetic_acid",
     "smiles": "OC(=O)Cc1c[nH]c2cc(F)ccc12",
     "rationale": "F→selectivity, indole→ARO, IAA→POS anchor"},

    # 3. ナフタレン系 (大きな疎水性界面向け)
    {"name": "2_naphthaleneacetic_acid",
     "smiles": "OC(=O)Cc1ccc2ccccc2c1",
     "rationale": "Naphthalene→large ARO/HYD pocket, COOH→POS"},
    {"name": "1_naphthylpropionic_acid",
     "smiles": "OC(=O)CCc1cccc2ccccc12",
     "rationale": "Extended naphthalene+acid for deep pocket"},

    # 4. スルホンアミド系 (COOHより強いLYS/ARG結合)
    {"name": "phenylsulfonamide",
     "smiles": "NS(=O)(=O)c1ccccc1",
     "rationale": "SO2NH2→LYS/ARG(强POS anchor), phenyl→ARO/HYD"},
    {"name": "benzenesulfonic_acid",
     "smiles": "OS(=O)(=O)c1ccccc1",
     "rationale": "SO3H→very strong POS anchor, phenyl→HYD"},

    # 5. テトラゾール系 (COOHのバイオアイソスター、膜透過性高い)
    {"name": "phenyl_tetrazole",
     "smiles": "c1ccc(-c2nnn[nH]2)cc1",
     "rationale": "Tetrazole(bioisostere of COOH)→LYS/ARG, phenyl→ARO"},
    {"name": "benzyl_tetrazole",
     "smiles": "c1ccc(Cc2nnn[nH]2)cc1",
     "rationale": "Tetrazole→POS anchor, benzyl→HYD+flexible"},

    # 6. ビフェニル・拡張芳香族系
    {"name": "biphenyl_4_carboxylic_acid",
     "smiles": "OC(=O)c1ccc(-c2ccccc2)cc1",
     "rationale": "Two phenyl rings→TYR pi-stack+HYD, COOH→POS"},
    {"name": "4_hydroxyphenylacetic_acid",
     "smiles": "OC(=O)Cc1ccc(O)cc1",
     "rationale": "OH→TYR H-bond, phenyl→ARO, COOH→LYS/ARG"},

    # 7. 複合機能基 (POS + HYD + POL)
    {"name": "tryptophan_analog",
     "smiles": "NC(Cc1c[nH]c2ccccc12)C(=O)O",
     "rationale": "Indole→ARO/HYD, NH2→NEG pocket, COOH→POS anchor"},
    {"name": "phenylalanine_analog",
     "smiles": "NC(Cc1ccccc1)C(=O)O",
     "rationale": "Phenyl→ARO/HYD, NH2+COOH zwitterion→both POS/NEG"},
    {"name": "methylindole_propanoic",
     "smiles": "OC(=O)CCc1c[nH]c2cc(C)ccc12",
     "rationale": "Methyl-indole→larger HYD, propanoate→POS anchor"},

    # 8. BRICS組み立て系 (フラグメント連結)
    {"name": "cyclohexyl_benzoic",
     "smiles": "OC(=O)c1ccc(C2CCCCC2)cc1",
     "rationale": "Cyclohexyl→LEU/ILE HYD, benzoate→POS anchor"},
    {"name": "fluorobenzyl_tetrazole",
     "smiles": "c1cc(F)ccc1Cc1nnn[nH]1",
     "rationale": "F-benzyl→selective HYD, tetrazole→POS anchor"},
    {"name": "naphthalene_sulfonamide",
     "smiles": "NS(=O)(=O)c1ccc2ccccc2c1",
     "rationale": "Naphthalene→large ARO/HYD interface, SO2NH2→POS"},

    # 9. 多価アニオン (複数POS残基に対応)
    {"name": "phthalic_acid",
     "smiles": "OC(=O)c1ccccc1C(=O)O",
     "rationale": "Two COOH→LYS48+LYS46 dual anchor, phenyl→HYD"},
    {"name": "isophthalic_acid",
     "smiles": "OC(=O)c1cccc(C(=O)O)c1",
     "rationale": "1,3-diCOOH→two POS residues, phenyl core→HYD"},
    {"name": "fumaric_acid_phenyl",
     "smiles": "OC(=O)/C=C/c1ccccc1",
     "rationale": "Trans-alkene→planar, phenyl→ARO, COOH→POS"},
]


# ──────────────────────────────────────────
# ポケット残基スコアリング
# ──────────────────────────────────────────

def score_protein_pocket(pdb_path: str,
                          contacts_json: str,
                          sasa_csv: str,
                          chain_protein: str = "A",
                          chain_peptide: str = "B") -> dict:
    """
    3手法の統合スコアでタンパク質側のポケット残基をランキングする。

    スコア構成 (正規化して合算):
      - 距離ベース接触数 × 0.4
      - ΔSASA (タンパク質側) × 0.4
      - PRODIGY Cβ接触 × 0.2

    Returns:
        {(resnum, resname): {"score": float, "type": str,
                              "n_contacts": int, "delta_sasa": float}}
    """
    # 1. 距離ベース接触数
    with open(contacts_json) as f:
        contacts = json.load(f)

    dist_counts = defaultdict(int)
    for c in contacts:
        key = (c["protein_resnum"], c["protein_res"])
        dist_counts[key] += 1

    # 2. ΔSASA (タンパク質側)
    sasa_vals = {}
    with open(sasa_csv) as f:
        for row in csv.DictReader(f):
            if row["chain"] == chain_protein:
                key = (int(row["resnum"]), row["resname"])
                sasa_vals[key] = float(row["delta_sasa"])

    # 3. PRODIGY Cβ接触 (タンパク質側の関与残基)
    prodigy_counts = defaultdict(int)
    parser = PDBParser(QUIET=True)
    struct  = parser.get_structure("cplx", pdb_path)
    model   = struct[0]

    def _get_cb(res):
        for a in res.get_atoms():
            if a.get_name() == "CB": return a
        for a in res.get_atoms():
            if a.get_name() == "CA": return a
        return None

    prot_res = list(model[chain_protein].get_residues())
    pep_res  = list(model[chain_peptide].get_residues())
    for pr in prot_res:
        cb_p = _get_cb(pr)
        if cb_p is None: continue
        for pe in pep_res:
            cb_e = _get_cb(pe)
            if cb_e is None: continue
            if float(np.linalg.norm(cb_p.coord - cb_e.coord)) <= 5.5:
                key = (pr.get_id()[1], pr.get_resname())
                prodigy_counts[key] += 1

    # 全キーを統合
    all_keys = set(dist_counts) | set(sasa_vals) | set(prodigy_counts)

    # 正規化
    max_d = max(dist_counts.values(), default=1)
    max_s = max(sasa_vals.values(),   default=1)
    max_p = max(prodigy_counts.values(), default=1)

    pocket = {}
    for key in all_keys:
        rnum, rname = key
        d = dist_counts.get(key, 0) / max_d
        s = sasa_vals.get(key, 0)   / max_s
        p = prodigy_counts.get(key, 0) / max_p
        score = 0.4 * d + 0.4 * s + 0.2 * p
        pocket[key] = {
            "score"     : round(score, 4),
            "type"      : RES_TYPE.get(rname, "UNK"),
            "n_contacts": dist_counts.get(key, 0),
            "delta_sasa": round(sasa_vals.get(key, 0), 2),
            "prodigy_n" : prodigy_counts.get(key, 0),
        }

    return dict(sorted(pocket.items(), key=lambda x: -x[1]["score"]))


def get_pocket_coords(pdb_path: str, pocket_residues: list,
                      chain_protein: str = "A") -> np.ndarray:
    """ポケット残基の重心座標を計算する"""
    parser = PDBParser(QUIET=True)
    struct = parser.get_structure("cplx", pdb_path)
    model  = struct[0]

    pocket_res_nums = {r for r, _ in pocket_residues}
    coords = []
    for res in model[chain_protein].get_residues():
        if res.get_id()[1] in pocket_res_nums:
            for atom in res.get_atoms():
                if atom.element != "H":
                    coords.append(atom.coord)
    return np.array(coords) if coords else np.zeros((1, 3))


# ──────────────────────────────────────────
# 相補的ファーマコフォア導出
# ──────────────────────────────────────────

def derive_complementary_pharmacophore(pocket: dict,
                                        top_n: int = 8) -> dict:
    """
    上位N残基のタイプ分布から相補的なリガンド特徴を導出する。

    Returns:
        {"required": [feat_type, ...], "type_counts": Counter,
         "dominant_type": str, "fragments_needed": dict}
    """
    top_residues = list(pocket.items())[:top_n]
    type_counts  = Counter(v["type"] for _, v in top_residues)
    type_scores  = defaultdict(float)
    for (_, name), v in top_residues:
        type_scores[v["type"]] += v["score"]

    # 相補的リガンド特徴
    complementary_map = {
        "POS": "NEG_anchor",    # 正電荷残基 → 負電荷アンカー
        "NEG": "POS_anchor",    # 負電荷残基 → 正電荷アンカー
        "ARO": "ARO_scaffold",  # 芳香族 → 芳香族スタッキング
        "HYD": "HYD_scaffold",  # 疎水性 → 疎水性コア
        "POL": "HBOND_group",   # 極性 → H結合基
    }

    required = []
    # スコア加重でソート → 最重要な相補的特徴を決定
    for rtype, _ in sorted(type_scores.items(), key=lambda x: -x[1]):
        feat = complementary_map.get(rtype)
        if feat and feat not in required:
            required.append(feat)

    dominant = max(type_scores, key=type_scores.get, default="HYD")

    return {
        "required"        : required,
        "type_counts"     : type_counts,
        "type_scores"     : dict(type_scores),
        "dominant_type"   : dominant,
        "complementary"   : complementary_map.get(dominant, "HYD_scaffold"),
    }


# ──────────────────────────────────────────
# 分子生成
# ──────────────────────────────────────────

def filter_by_pharmacophore(mol, required: list) -> bool:
    """
    分子が必要なファーマコフォア特徴を持つかチェック。
    """
    smarts_map = {
        "NEG_anchor" : ["[CX3](=O)[OH]",       # COOH
                        "[SX4](=O)(=O)[NH2]",   # SO2NH2
                        "c1nnn[nH]1",            # tetrazole
                        "[PX4](=O)([OH])[OH]"], # phosphonate
        "POS_anchor" : ["[NH2]", "[NH3+]", "NC(=N)N"],
        "ARO_scaffold": ["c1ccccc1", "c1ccncc1", "c1ccc2[nH]ccc2c1"],
        "HYD_scaffold": ["[CX4H2][CX4H2]", "C1CCCCC1", "c1ccccc1"],
        "HBOND_group" : ["[OH]", "[NH2]", "C(=O)N"],
    }
    for feat in required:
        patterns = smarts_map.get(feat, [])
        if not patterns:
            continue
        matched = any(mol.HasSubstructMatch(Chem.MolFromSmarts(p))
                      for p in patterns if Chem.MolFromSmarts(p))
        if not matched:
            return False
    return True


def generate_brics_molecules(seed_smiles: list, n_mols: int = 20) -> list:
    """
    BRICSを使って種分子からフラグメントを分解し、再組み立てして新分子を生成。
    """
    from rdkit.Chem import BRICS

    fragments = set()
    for smi in seed_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        frags = BRICS.BRICSDecompose(mol)
        fragments.update(frags)

    if not fragments:
        return []

    frag_mols = [Chem.MolFromSmiles(f) for f in fragments]
    frag_mols = [m for m in frag_mols if m is not None]

    try:
        new_mols = list(BRICS.BRICSBuild(frag_mols))[:n_mols]
        return [Chem.MolToSmiles(m) for m in new_mols if m is not None]
    except Exception:
        return []


def calculate_props(mol) -> dict:
    """分子プロパティ計算 (utils.drug_likeness に委譲)"""
    from utils.drug_likeness import calculate_drug_likeness
    return calculate_drug_likeness(mol)


def generate_3d(mol) -> Optional[object]:
    """3D配座生成"""
    mol_h = Chem.AddHs(mol)
    params = AllChem.EmbedParameters()
    params.randomSeed = 42
    if AllChem.EmbedMolecule(mol_h, params) < 0:
        # 失敗したら ETKDG で再試行
        AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
    if mol_h.GetNumConformers() == 0:
        return None
    AllChem.MMFFOptimizeMolecule(mol_h)
    return Chem.RemoveHs(mol_h)


def assemble_candidates(pocket: dict,
                         pharmacophore: dict,
                         top_pocket_n: int = 8) -> list[dict]:
    """
    ポケット情報に基づいて候補分子を組み立てる。
    1. テンプレートライブラリから適合分子を選択
    2. BRICSで新規分子を生成
    3. フィルタリング (Ro5 + ファーマコフォア適合)
    """
    required = pharmacophore["required"]
    candidates = []

    # ── 1. テンプレートから選択 ──
    print("  [1] テンプレートライブラリから選択...")
    for tmpl in MOLECULE_TEMPLATES:
        mol = Chem.MolFromSmiles(tmpl["smiles"])
        if mol is None:
            continue
        props = calculate_props(mol)
        if not props["Ro5"]:
            continue
        pharm_ok = filter_by_pharmacophore(mol, required[:2])  # 必須特徴の最初の2つ
        candidates.append({
            "name"      : tmpl["name"],
            "smiles"    : Chem.MolToSmiles(mol),
            "mol"       : mol,
            "rationale" : tmpl["rationale"],
            "source"    : "template",
            "pharm_ok"  : pharm_ok,
            **props,
        })

    # ── 2. BRICS組み立て ──
    print("  [2] BRICS フラグメント組み立て...")
    seed_smiles = [t["smiles"] for t in MOLECULE_TEMPLATES[:8]]
    brics_smiles = generate_brics_molecules(seed_smiles, n_mols=30)
    for smi in brics_smiles:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        props = calculate_props(mol)
        if not props["Ro5"]:
            continue
        if props["MW"] < 150 or props["MW"] > 450:
            continue
        pharm_ok = filter_by_pharmacophore(mol, required[:1])
        name = f"brics_{len([c for c in candidates if c['source']=='brics'])}"
        candidates.append({
            "name"      : name,
            "smiles"    : smi,
            "mol"       : mol,
            "rationale" : "BRICS assembly from template fragments",
            "source"    : "brics",
            "pharm_ok"  : pharm_ok,
            **props,
        })

    # ── 3. フィルタリングと重複除去 ──
    seen   = set()
    unique = []
    for c in candidates:
        canonical = Chem.MolToSmiles(c["mol"])
        if canonical in seen:
            continue
        seen.add(canonical)
        unique.append({**c, "smiles": canonical})

    # ファーマコフォア適合のものを優先
    unique.sort(key=lambda x: (not x["pharm_ok"], -x["MW"]))
    print(f"  → {len(unique)} 候補分子生成 "
          f"(内ファーマコフォア適合: {sum(c['pharm_ok'] for c in unique)})")
    return unique


# ──────────────────────────────────────────
# 表示・保存
# ──────────────────────────────────────────

def print_pocket_summary(pocket: dict, pharmacophore: dict, top_n: int = 10):
    print("\n" + "=" * 70)
    print("  タンパク質結合ポケット 統合解析")
    print("=" * 70)
    print(f"  {'残基':<10} {'スコア':>6}  {'タイプ':<5} "
          f"{'接触数':>6}  {'ΔSASA(Å²)':>10}  {'Cβ接触':>6}")
    print("  " + "-" * 55)
    for (rnum, rname), v in list(pocket.items())[:top_n]:
        bar = "█" * int(v["score"] * 20)
        print(f"  {rname}{rnum:<7} {v['score']:>6.3f}  {v['type']:<5} "
              f"{v['n_contacts']:>6}  {v['delta_sasa']:>10.1f}  "
              f"{v['prodigy_n']:>6}  {bar}")
    print()
    print("  相補的ファーマコフォア要件:")
    for feat in pharmacophore["required"]:
        print(f"    ▶ {feat}")
    print(f"  ポケット支配タイプ: {pharmacophore['dominant_type']}")
    print("=" * 70)


def save_candidates_sdf(candidates: list[dict], output_path: str):
    writer = Chem.SDWriter(output_path)
    saved = 0
    for c in candidates:
        mol_3d = generate_3d(c["mol"])
        target = mol_3d if mol_3d else c["mol"]
        if target is None:
            continue
        target.SetProp("_Name",    c["name"])
        target.SetProp("SMILES",   c["smiles"])
        target.SetProp("Source",   c["source"])
        target.SetProp("Rationale",c["rationale"])
        target.SetProp("PharmOK",  str(c["pharm_ok"]))
        for k in ("MW", "LogP", "HBD", "HBA", "PSA", "RotBonds"):
            target.SetProp(k, str(c.get(k, "")))
        writer.write(target)
        saved += 1
    writer.close()
    print(f"  候補分子 SDF 保存: {output_path} ({saved}分子)")


def save_pocket_csv(pocket: dict, output_path: str):
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "resnum", "resname", "type",
                         "score", "n_contacts", "delta_sasa", "prodigy_n"])
        for rank, ((rnum, rname), v) in enumerate(pocket.items(), 1):
            writer.writerow([rank, rnum, rname, v["type"],
                             v["score"], v["n_contacts"],
                             v["delta_sasa"], v["prodigy_n"]])
    print(f"  ポケット残基 CSV 保存: {output_path}")


def print_candidates_table(candidates: list[dict]):
    print(f"\n  {'名前':<30} {'MW':>6} {'LogP':>6} {'PSA':>6} "
          f"{'Ro5'} {'ファーマコフォア'} {'出典'}")
    print("  " + "-" * 80)
    for c in candidates[:25]:
        flag = "✓" if c["Ro5"]      else "✗"
        phok = "✓" if c["pharm_ok"] else "△"
        print(f"  {c['name']:<30} {c['MW']:>6.0f} {c['LogP']:>6.2f} "
              f"{c['PSA']:>6.1f} {flag}   {phok}          {c['source']}")


def plot_pocket_features(pocket: dict, pharmacophore: dict, output_path: str):
    """ポケット残基タイプ分布の可視化"""
    import plot_utils  # noqa: F401
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    TYPE_COLOR = {"POS": "#e74c3c", "NEG": "#3498db", "ARO": "#e67e22",
                  "HYD": "#f1c40f", "POL": "#2ecc71", "GLY": "#95a5a6",
                  "UNK": "#bdc3c7"}

    top = list(pocket.items())[:12]
    labels = [f"{name}{num}" for (num, name), _ in top]
    scores = [v["score"]  for _, v in top]
    colors = [TYPE_COLOR.get(v["type"], "#bdc3c7") for _, v in top]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # 棒グラフ
    bars = ax1.barh(labels[::-1], scores[::-1], color=colors[::-1],
                    edgecolor="white")
    ax1.set_xlabel("統合スコア", fontsize=11)
    ax1.set_title("タンパク質ポケット残基 統合スコア", fontsize=12, fontweight="bold")

    patches = [mpatches.Patch(color=c, label=l)
               for l, c in TYPE_COLOR.items() if l != "UNK"]
    ax1.legend(handles=patches, loc="lower right", fontsize=8)

    # タイプ分布円グラフ
    tc = pharmacophore["type_scores"]
    tc_filt = {k: v for k, v in tc.items() if v > 0}
    pie_colors = [TYPE_COLOR.get(k, "#bdc3c7") for k in tc_filt]
    ax2.pie(tc_filt.values(), labels=tc_filt.keys(),
            colors=pie_colors, autopct="%1.0f%%",
            textprops={"fontsize": 11})
    ax2.set_title("ポケット残基タイプ分布", fontsize=12, fontweight="bold")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  ポケット特性グラフ 保存: {output_path}")


# ──────────────────────────────────────────
# メイン
# ──────────────────────────────────────────

def run(pdb_path       : str = "Protein_Peptide.pdb",
        contacts_json  : str = "results/contacts.json",
        sasa_csv       : str = "results/sasa_analysis.csv",
        output_dir     : str = "results",
        top_pocket_n   : int = 8,
        chain_protein  : str = "A",
        chain_peptide  : str = "B") -> dict:

    os.makedirs(output_dir, exist_ok=True)

    print("\n[ポケット解析] タンパク質側接触残基を統合スコアリング...")
    pocket = score_protein_pocket(
        pdb_path, contacts_json, sasa_csv,
        chain_protein, chain_peptide
    )

    print("[ファーマコフォア] 相補的特徴を導出...")
    pharmacophore = derive_complementary_pharmacophore(pocket, top_n=top_pocket_n)

    print_pocket_summary(pocket, pharmacophore)
    save_pocket_csv(pocket, os.path.join(output_dir, "pocket_residues.csv"))
    plot_pocket_features(pocket, pharmacophore,
                         os.path.join(output_dir, "pocket_features.png"))

    print("\n[分子生成] 相補的低分子を生成中...")
    candidates = assemble_candidates(pocket, pharmacophore, top_pocket_n)
    print_candidates_table(candidates)

    sdf_path = os.path.join(output_dir, "candidate_pocket_ligands.sdf")
    save_candidates_sdf(candidates, sdf_path)

    return {
        "pocket"       : pocket,
        "pharmacophore": pharmacophore,
        "candidates"   : candidates,
        "sdf_path"     : sdf_path,
    }


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--pdb",          default="Protein_Peptide.pdb")
    p.add_argument("--contacts",     default="results/contacts.json")
    p.add_argument("--sasa",         default="results/sasa_analysis.csv")
    p.add_argument("--output-dir",   default="results")
    p.add_argument("--top-pocket-n", type=int, default=8)
    p.add_argument("--dock",         action="store_true")
    args = p.parse_args()

    result = run(
        pdb_path      = args.pdb,
        contacts_json = args.contacts,
        sasa_csv      = args.sasa,
        output_dir    = args.output_dir,
        top_pocket_n  = args.top_pocket_n,
    )

    if args.dock:
        print("\n[ドッキング] smina でドッキング実行...")
        from dock_with_smina import (extract_chain, run_smina_per_molecule,
                                     print_docking_summary, save_docking_results,
                                     plot_docking_scores, rescore_all_docked,
                                     save_rescore_csv, plot_score_decomposition)

        dock_dir = os.path.join(args.output_dir, "docking_pocket")
        os.makedirs(dock_dir, exist_ok=True)

        receptor_pdb = os.path.join(dock_dir, "receptor.pdb")
        peptide_ref  = os.path.join(dock_dir, "peptide_ref.pdb")
        extract_chain(args.pdb, "A", receptor_pdb)
        extract_chain(args.pdb, "B", peptide_ref)

        dock_results = run_smina_per_molecule(
            receptor_pdb, result["sdf_path"], peptide_ref, dock_dir,
            exhaustiveness=8, num_modes=9,
        )
        print_docking_summary(dock_results)
        save_docking_results(dock_results,
                             os.path.join(dock_dir, "docking_results.csv"))
        plot_docking_scores(dock_results,
                            os.path.join(dock_dir, "docking_scores.png"))

        print("\n[Rescoring] smina score_only でスコア分解...")
        enriched = rescore_all_docked(dock_results, receptor_pdb, dock_dir)
        save_rescore_csv(enriched, os.path.join(dock_dir, "rescore_terms.csv"))
        plot_score_decomposition(enriched,
                                 os.path.join(dock_dir, "score_decomposition.png"))
