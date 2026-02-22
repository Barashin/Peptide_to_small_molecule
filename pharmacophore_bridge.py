"""
pharmacophore_bridge.py
=======================
2つのタンパク質ポケット残基を指定し、そのCβ–Cβ距離を目標に
ファーマコフォア拘束付き3D低分子を生成するモジュール。

アルゴリズム概要:
  1. 2残基のCβ座標をPDBから取得し、目標Cβ–Cβ距離を計算
  2. 各残基タイプに相補的なアンカー基 (Anchor A, Anchor B) を選択
  3. 距離に合わせたリンカー長を選択
  4. SMILES組み立て:  AnchorA + Linker + AnchorB
  5. RDKit距離幾何学:
       GetMoleculeBoundsMatrix() → 距離制約を設定 →
       DoTriangleSmoothing()    → EmbedMolecule (制約付き3D埋め込み)
  6. MMFFで構造最適化
  7. SDF保存 → smina でドッキング

使用例:
  python pharmacophore_bridge.py --point1 LYS48 --point2 LEU50
  python pharmacophore_bridge.py --point1 LYS46 --point2 TYR49 --no-dock
  python pharmacophore_bridge.py --list-pairs

CLI引数:
  --point1 RESNAME+NUM  残基1 (例: LYS48)
  --point2 RESNAME+NUM  残基2 (例: LEU50)
  --chain   STR         タンパク質チェーンID (デフォルト: A)
  --no-dock             ドッキングをスキップ
  --n-confs INT         コンフォマー数 (デフォルト: 50)
  --list-pairs          既知ペアを一覧表示して終了
"""

import os
import re
import sys
import json
import argparse
import subprocess
import itertools
import csv
import numpy as np
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors
from rdkit.Chem import rdDistGeom

import plot_utils  # noqa: F401  (Agg バックエンド + 日本語フォント)
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

import linker_library as _ll

# FEgrow リンカー DB の遅延ロード
_LINKER_DB: list | None = None

def _get_linker_db() -> list:
    global _LINKER_DB
    if _LINKER_DB is None:
        _LINKER_DB = _ll.load_linker_db()
    return _LINKER_DB

# ──────────────────────────────────────────
# パス設定
# ──────────────────────────────────────────
BASE_DIR       = Path(__file__).parent
PDB_PATH       = str(BASE_DIR / "Protein_Peptide.pdb")
SMINA_PATH     = str(BASE_DIR / "smina.osx.12")
RECEPTOR_PATH  = str(BASE_DIR / "results" / "docking" / "receptor.pdb")
PEPTIDE_REF    = str(BASE_DIR / "results" / "docking" / "peptide_ref.pdb")
OUT_DIR        = BASE_DIR / "results" / "bridge"

PROTEIN_CHAIN  = "A"

# ──────────────────────────────────────────
# 残基タイプ定義
# ──────────────────────────────────────────
from utils.residue_defs import RES_TYPE

# ──────────────────────────────────────────
# アンカー基ライブラリ
#
# 各エントリ:
#   name         : アンカー名
#   smiles       : SMILES文字列 (A側は左端がファーマコフォア,
#                                B側は左端がリンカーとの接続点)
#   pharm_offset : ファーマコフォア原子のSMILES内重原子インデックス
#   desc         : 説明
#
# A側 (prefix): smiles + linker_smiles + B_smiles
# B側 (suffix): A_smiles + linker_smiles + smiles
# ──────────────────────────────────────────

# A側アンカー: 「タンパク質残基と相互作用する基」が左端に来るSMILES
ANCHOR_A = {
    "POS": [
        {
            "name": "carboxylic_acid",
            "smiles": "OC(=O)",   # O(0)–C(1)(=O(2))–[linker]
            "pharm_offset": 0,    # OH酸素 (LYS ε-NH3+と塩橋)
            "desc": "カルボン酸 (LYS/ARG NH₃⁺と塩橋)",
        },
        {
            "name": "sulfonamide",
            "smiles": "NS(=O)(=O)",  # N(0)–S(1)(=O(2))(=O(3))–[linker]
            "pharm_offset": 0,
            "desc": "スルホンアミド (LYS/ARG NH₃⁺と水素結合)",
        },
        {
            "name": "tetrazole",
            "smiles": "c1nnn[nH]1",  # テトラゾール 5原子, 最後のC(4)がlinker側
            "pharm_offset": 0,
            "desc": "テトラゾール (COOH等価体, 膜透過性高)",
        },
    ],
    "NEG": [
        {
            "name": "amine",
            "smiles": "N",
            "pharm_offset": 0,
            "desc": "一級アミン (ASP/GLU COO⁻と塩橋)",
        },
        {
            "name": "guanidine",
            "smiles": "NC(=N)N",  # N(0)–C(1)(=N(2))–N(3)–[linker]
            "pharm_offset": 0,
            "desc": "グアニジン (ASP/GLU COO⁻と強塩橋)",
        },
    ],
    "ARO": [
        {
            "name": "phenyl",
            "smiles": "c1ccccc1",  # 6原子
            "pharm_offset": 0,
            "desc": "フェニル (π–πスタッキング)",
        },
    ],
    "HYD": [
        {
            "name": "isobutyl",
            "smiles": "CC(C)C",   # C(0)–C(1)(C(2))–C(3)–[linker]
            "pharm_offset": 1,    # 分岐C
            "desc": "イソブチル (LEU/ILE疎水接触)",
        },
        {
            "name": "cyclohexylmethyl",
            "smiles": "C1CCCCC1C",  # 7原子, 最後のC(6)がlinker側
            "pharm_offset": 0,
            "desc": "シクロヘキシルメチル (LEU/ILE疎水コア)",
        },
    ],
    "POL": [
        {
            "name": "hydroxymethyl",
            "smiles": "OC",    # O(0)–C(1)–[linker]
            "pharm_offset": 0,
            "desc": "ヒドロキシメチル (HBA/HBD)",
        },
        {
            "name": "amide",
            "smiles": "NC(=O)",  # N(0)–C(1)(=O(2))–[linker]
            "pharm_offset": 0,
            "desc": "アミド",
        },
    ],
}

# B側アンカー: 「リンカー接続点」が左端, 「ファーマコフォア」が右側
ANCHOR_B = {
    "POS": [
        {
            "name": "carboxylic_acid",
            "smiles": "C(=O)O",   # [linker]–C(0)(=O(1))–O(2)H
            "pharm_offset": 2,    # OH酸素
            "desc": "カルボン酸 (LYS/ARG NH₃⁺と塩橋)",
        },
    ],
    "NEG": [
        {
            "name": "amine",
            "smiles": "N",
            "pharm_offset": 0,
            "desc": "アミン (ASP/GLU COO⁻と塩橋)",
        },
    ],
    "ARO": [
        {
            "name": "phenyl",
            "smiles": "c1ccccc1",
            "pharm_offset": 0,
            "desc": "フェニル (π–πスタッキング)",
        },
        {
            "name": "naphthalene",
            "smiles": "c1ccc2ccccc2c1",  # 10原子
            "pharm_offset": 0,
            "desc": "ナフタレン (拡張π–πスタッキング)",
        },
        {
            "name": "indole",
            "smiles": "c1ccc2[nH]ccc2c1",  # 9原子
            "pharm_offset": 0,
            "desc": "インドール (TRP mimetic)",
        },
    ],
    "HYD": [
        {
            "name": "cyclohexyl",
            "smiles": "C1CCCCC1",  # 6原子
            "pharm_offset": 0,
            "desc": "シクロヘキシル (LEU/ILE疎水コア)",
        },
        {
            "name": "phenyl",
            "smiles": "c1ccccc1",
            "pharm_offset": 0,
            "desc": "フェニル (疎水性/π)",
        },
        {
            "name": "isopropyl",
            "smiles": "C(C)C",   # [linker]–C(0)(–C(1))–C(2)
            "pharm_offset": 1,
            "desc": "イソプロピル (ALA/VAL疎水接触)",
        },
    ],
    "POL": [
        {
            "name": "hydroxyl",
            "smiles": "O",
            "pharm_offset": 0,
            "desc": "ヒドロキシル (HBA)",
        },
        {
            "name": "amide",
            "smiles": "C(=O)N",   # [linker]–C(0)(=O(1))–N(2)
            "pharm_offset": 2,
            "desc": "アミド N (HBD)",
        },
    ],
}

# ──────────────────────────────────────────
# リンカーライブラリは linker_library.py (FEgrow) に委譲
# ──────────────────────────────────────────

# ──────────────────────────────────────────
# 既知ポケットペア (解析済み, 参考用)
# ──────────────────────────────────────────
KNOWN_PAIRS = [
    {"res1": "LYS48", "res2": "LEU50", "dist": 7.3, "type": "POS–HYD"},
    {"res1": "LYS46", "res2": "TYR49", "dist": 6.6, "type": "POS–ARO"},
    {"res1": "ILE21", "res2": "TYR25", "dist": 6.1, "type": "HYD–ARO"},
    {"res1": "ARG28", "res2": "TYR49", "dist": 9.1, "type": "POS–ARO"},
    {"res1": "LYS46", "res2": "LEU50", "dist": 10.0,"type": "POS–HYD"},
]


# ══════════════════════════════════════════
# ユーティリティ関数
# ══════════════════════════════════════════

def parse_residue_id(res_str: str) -> tuple[str, int]:
    """'LYS48' → ('LYS', 48)"""
    m = re.match(r"([A-Za-z]+)(\d+)", res_str.strip())
    if not m:
        raise ValueError(f"残基ID解析エラー: '{res_str}' (例: LYS48)")
    return m.group(1).upper(), int(m.group(2))


def get_cbeta_coords(pdb_path: str, chain_id: str, resnum: int) -> tuple[np.ndarray, str]:
    """指定残基のCβ (GLY→Cα) 座標とresname を返す"""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                if residue.id[1] == resnum:
                    resname = residue.resname.strip()
                    for atom_name in ("CB", "CA"):
                        if atom_name in residue:
                            return np.array(residue[atom_name].get_coord()), resname
                    raise ValueError(f"Cβ/Cα が残基 {chain_id}:{resnum} に見つかりません")
    raise ValueError(f"残基 {chain_id}:{resnum} が {pdb_path} に存在しません")


def get_n_atoms(smiles: str) -> int:
    """SMILES の重原子数を返す (空文字列→0)"""
    if not smiles:
        return 0
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumHeavyAtoms() if mol else 0


def select_linkers(target_dist: float, tolerance: float = 2.5) -> list[dict]:
    """
    target_dist ± tolerance (Å) に収まる FEgrow リンカーを返す。
    該当なしの場合は最近傍のリンカーを1つ返す (linker_library が保証)。
    """
    return _ll.select_linkers_by_dist(
        _get_linker_db(), target_dist, tolerance=tolerance, max_n=5,
    )


def build_bridge_smiles(
    anchor_a: dict,
    linker: dict,
    anchor_b: dict,
) -> tuple[str, int, int] | None:
    """
    FEgrow リンカーを使ってブリッジ SMILES を構築し、ファーマコフォア原子インデックスを返す。

    アンカー官能基 (carboxylic acid, phenyl 等) が異種原子で始まる場合の原子価エラーを
    回避するため、線形コア抽出 + 文字列連結方式 (build_bridge_smiles_linear) を使用する。
    非線形リンカーは None を返し、generate_bridges() でスキップされる。

    Returns:
        (smiles, pharm_a_idx, pharm_b_idx) または None (無効SMILESまたは非線形リンカー)
    """
    return _ll.build_bridge_smiles_linear(anchor_a, linker, anchor_b)


def embed_with_constraint(
    smiles: str,
    pharm_a_idx: int,
    pharm_b_idx: int,
    target_dist: float,
    tolerance: float = 2.5,
    n_confs: int = 50,
    seed: int = 42,
) -> tuple[Chem.Mol | None, list[int]]:
    """
    RDKit 距離幾何学で制約付きコンフォマーを生成する。

    手順:
      1. SMILES → Mol → AddHs
      2. GetMoleculeBoundsMatrix() でデフォルトバウンズ行列取得
      3. pharm_a ↔ pharm_b 間の上下限を target_dist ± tolerance に設定
      4. DoTriangleSmoothing() で三角不等式を満足させる
      5. EmbedMolecule() で制約付き3D埋め込み
      6. MMFF最適化

    Returns:
        (mol_with_confs, list_of_conf_ids)
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, []

    mol = Chem.RWMol(mol)
    AllChem.Compute2DCoords(mol)
    mol_h = Chem.AddHs(mol)

    # ── バウンズ行列の取得と制約設定 ──
    bm = rdDistGeom.GetMoleculeBoundsMatrix(mol_h)
    i, j = sorted([pharm_a_idx, pharm_b_idx])

    # 上三角 = 上限距離, 下三角 = 下限距離
    upper = min(target_dist + tolerance, bm[i][j])
    lower = max(target_dist - tolerance, bm[j][i], 1.0)
    if lower < upper:
        bm[i][j] = upper
        bm[j][i] = lower
        use_bm = True
    else:
        print("    ⚠ 距離制約の上下限が矛盾。拘束なしで埋め込みます。")
        use_bm = False

    # ── コンフォマー生成 ──
    conf_ids: list[int] = []

    if use_bm:
        # EmbedParameters.SetBoundsMat (rdkit 2022+)
        try:
            params = rdDistGeom.EmbedParameters()
            params.randomSeed            = seed
            params.useRandomCoords       = True
            params.ignoreSmoothingFailures = True
            params.SetBoundsMat(bm)
            conf_ids = list(rdDistGeom.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params))
        except Exception as e:
            print(f"    ⚠ SetBoundsMat 失敗 ({e})。拘束なしで再試行。")
            use_bm = False

    # 拘束なしフォールバック
    if not conf_ids:
        print("    ⚠ 拘束なしで埋め込みます。")
        params2 = rdDistGeom.EmbedParameters()
        params2.randomSeed      = seed
        params2.useRandomCoords = True
        conf_ids = list(rdDistGeom.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params2))

    if not conf_ids:
        conf_ids = list(AllChem.EmbedMultipleConfs(
            mol_h, numConfs=n_confs, randomSeed=seed, useRandomCoords=True,
        ))

    if not conf_ids:
        return None, []

    # ── MMFF 最適化 ──
    try:
        results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, numThreads=0, maxIters=500)
    except Exception:
        results = [(1, 0.0)] * len(conf_ids)

    return mol_h, conf_ids


def select_best_conformer(
    mol_h: Chem.Mol,
    conf_ids: list[int],
    pharm_a_idx: int,
    pharm_b_idx: int,
    target_dist: float,
) -> int:
    """
    コンフォマーの中から最良のものを選ぶ。
    基準: |実測距離 - 目標距離|² + MMFF エネルギー (正規化)

    Returns:
        最良コンフォマーの conf_id
    """
    if len(conf_ids) == 1:
        return conf_ids[0]

    scores = []
    for cid in conf_ids:
        conf = mol_h.GetConformer(cid)
        pos_a = np.array(conf.GetAtomPosition(pharm_a_idx))
        pos_b = np.array(conf.GetAtomPosition(pharm_b_idx))
        actual_dist = float(np.linalg.norm(pos_a - pos_b))
        dist_penalty = (actual_dist - target_dist) ** 2

        # MMFFエネルギー取得
        try:
            ff = AllChem.MMFFGetMoleculeForceField(
                mol_h, AllChem.MMFFGetMoleculeProperties(mol_h), confId=cid
            )
            energy = ff.CalcEnergy() if ff else 0.0
        except Exception:
            energy = 0.0

        scores.append((cid, dist_penalty, energy))

    # ペナルティとエネルギーを正規化して合算
    penalties = [s[1] for s in scores]
    energies  = [s[2] for s in scores]
    p_range = max(penalties) - min(penalties) or 1.0
    e_range = max(energies)  - min(energies)  or 1.0

    best_cid = min(
        scores,
        key=lambda s: 0.7 * (s[1] / p_range) + 0.3 * (s[2] / e_range),
    )[0]
    return best_cid


def calc_drug_likeness(mol: Chem.Mol) -> dict:
    """Lipinski Ro5 + Veber ルールを計算 (utils.drug_likeness に委譲)"""
    from utils.drug_likeness import calculate_drug_likeness
    return calculate_drug_likeness(mol)


def save_sdf(
    mol_h: Chem.Mol,
    conf_id: int,
    name: str,
    props: dict,
    out_path: str,
) -> None:
    """コンフォマー1つを SDF に保存"""
    mol_no_h = Chem.RemoveHs(mol_h)
    # Hなし分子に最良コンフォマーの座標をコピー
    try:
        embed_mol = Chem.RWMol(mol_no_h)
        AllChem.EmbedMolecule(embed_mol, randomSeed=42)
        # 重原子座標だけコピー
        conf_h = mol_h.GetConformer(conf_id)
        conf_no_h = embed_mol.GetConformer()
        for atom in mol_no_h.GetAtoms():
            idx = atom.GetIdx()
            pos = conf_h.GetAtomPosition(idx)
            conf_no_h.SetAtomPosition(idx, pos)
        out_mol = embed_mol
    except Exception:
        out_mol = Chem.RWMol(mol_h)

    out_mol.SetProp("_Name", name)
    out_mol.SetProp("SMILES", Chem.MolToSmiles(mol_no_h))
    for k, v in props.items():
        out_mol.SetProp(str(k), str(v))

    writer = Chem.SDWriter(out_path)
    writer.write(out_mol)
    writer.close()


def run_docking(sdf_path: str, receptor_path: str, ref_path: str,
                out_dir: str, name: str, exhaustiveness: int = 8,
                num_modes: int = 9) -> tuple[float | None, str | None]:
    """smina でドッキングを実行し、スコアとdocked SDFパスを返す"""
    if not os.path.isfile(SMINA_PATH):
        print(f"  ⚠ smina が見つかりません: {SMINA_PATH}")
        return None, None

    out_sdf = os.path.join(out_dir, f"docked_{name}.sdf")
    log_txt  = os.path.join(out_dir, f"log_{name}.txt")

    cmd = [
        SMINA_PATH,
        "--receptor",       receptor_path,
        "--ligand",         sdf_path,
        "--autobox_ligand", ref_path,
        "--autobox_add",    "4",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes",      str(num_modes),
        "--out",            out_sdf,
        "--log",            log_txt,
        "--cpu",            "4",
        "--quiet",
    ]
    try:
        os.chmod(SMINA_PATH, 0o755)
        subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    except Exception as e:
        print(f"    ドッキングエラー: {e}")
        return None, None

    # ログからスコアを抽出
    score = None
    if os.path.isfile(log_txt):
        with open(log_txt) as f:
            for line in f:
                m = re.match(r"^\s*\d+\s+([-\d.]+)\s+[\d.]+\s+[\d.]+", line)
                if m:
                    score = float(m.group(1))
                    break

    if score is None or not os.path.isfile(out_sdf):
        return None, None
    return score, out_sdf


# ══════════════════════════════════════════
# メイン生成関数
# ══════════════════════════════════════════

def generate_bridges(
    res1_str: str,
    res2_str: str,
    pdb_path: str  = PDB_PATH,
    chain_id: str  = PROTEIN_CHAIN,
    dock: bool     = True,
    n_confs: int   = 50,
    tolerance: float = 2.5,
    exhaustiveness: int = 8,
) -> list[dict]:
    """
    2残基間を橋渡しする低分子を生成する。

    Args:
        res1_str : 残基1 (例: "LYS48")
        res2_str : 残基2 (例: "LEU50")
        pdb_path : タンパク質PDBファイルパス
        chain_id : タンパク質チェーンID
        dock     : smina でドッキングするか
        n_confs  : コンフォマー生成数

    Returns:
        results: 各候補分子の情報を含む dict のリスト
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── ① 残基情報・座標取得 ──
    resname1, resnum1 = parse_residue_id(res1_str)
    resname2, resnum2 = parse_residue_id(res2_str)

    print(f"\n{'='*60}")
    print(f"ファーマコフォアブリッジ生成")
    print(f"  残基1: {resname1}{resnum1}  残基2: {resname2}{resnum2}")

    coords1, rn1 = get_cbeta_coords(pdb_path, chain_id, resnum1)
    coords2, rn2 = get_cbeta_coords(pdb_path, chain_id, resnum2)
    target_dist = float(np.linalg.norm(coords1 - coords2))
    print(f"  Cβ–Cβ 距離: {target_dist:.2f} Å")

    type1 = RES_TYPE.get(rn1, "HYD")
    type2 = RES_TYPE.get(rn2, "HYD")
    print(f"  残基タイプ: {rn1}({type1}) — {rn2}({type2})")

    # ── ② アンカーとリンカーを選択 ──
    anchors_a = ANCHOR_A.get(type1, ANCHOR_A["HYD"])
    anchors_b = ANCHOR_B.get(type2, ANCHOR_B["HYD"])

    print(f"  アンカーA ({type1}): {[a['name'] for a in anchors_a]}")
    print(f"  アンカーB ({type2}): {[a['name'] for a in anchors_b]}")

    # アンカーサイズを考慮してリンカー目標距離をペアごとに計算
    combos: list[tuple] = []
    linker_names_seen: set[str] = set()
    for anchor_a in anchors_a:
        for anchor_b in anchors_b:
            a_reach = _ll.estimate_anchor_reach(anchor_a["smiles"])
            b_reach = _ll.estimate_anchor_reach(anchor_b["smiles"])
            linker_target = max(0.0, target_dist - a_reach - b_reach)
            linkers = select_linkers(linker_target, tolerance)
            for lk in linkers:
                combos.append((anchor_a, lk, anchor_b))
                linker_names_seen.add(lk["name"])

    print(f"  リンカー候補 (FEgrow): {sorted(linker_names_seen)}")
    print(f"  組み合わせ数: {len(combos)}")
    print(f"{'='*60}")

    results  = []
    seen_smi = set()
    mol_idx  = 0

    for anchor_a, linker, anchor_b in combos:
        built = build_bridge_smiles(anchor_a, linker, anchor_b)
        if built is None:
            continue
        smiles, pharm_a_idx, pharm_b_idx = built

        # SMILES の正規化と重複排除
        mol_tmp = Chem.MolFromSmiles(smiles)
        if mol_tmp is None:
            continue
        canon = Chem.MolToSmiles(mol_tmp)
        if canon in seen_smi:
            continue
        seen_smi.add(canon)

        # Ro5 フィルタ (MWのみ事前確認)
        mw = Descriptors.ExactMolWt(mol_tmp)
        if mw > 600:
            continue

        mol_idx += 1
        name = (f"bridge_{res1_str}_{res2_str}"
                f"_{anchor_a['name']}_{linker['name']}_{anchor_b['name']}")
        name = name.replace("(", "").replace(")", "").replace(" ", "_")[:60]

        print(f"\n[{mol_idx:03d}] {name}")
        print(f"  SMILES: {canon}")
        print(f"  MW={mw:.1f}, pharm_A={pharm_a_idx}, pharm_B={pharm_b_idx}")

        # ── ③ 制約付き3D埋め込み ──
        mol_h, conf_ids = embed_with_constraint(
            canon, pharm_a_idx, pharm_b_idx,
            target_dist, tolerance, n_confs, seed=42 + mol_idx,
        )
        if mol_h is None or not conf_ids:
            print(f"  ✗ コンフォマー生成失敗")
            continue

        # ── ④ 最良コンフォマー選択と距離確認 ──
        best_cid = select_best_conformer(
            mol_h, conf_ids, pharm_a_idx, pharm_b_idx, target_dist
        )
        conf = mol_h.GetConformer(best_cid)
        pos_a = np.array(conf.GetAtomPosition(pharm_a_idx))
        pos_b = np.array(conf.GetAtomPosition(pharm_b_idx))
        actual_dist = float(np.linalg.norm(pos_a - pos_b))
        print(f"  実測ファーマコフォア間距離: {actual_dist:.2f} Å (目標: {target_dist:.2f} Å)")

        # ── ④b コンフォメーション適合性スコア ──
        try:
            from rigid_scaffold_design import score_pharmacophore_conformance, count_rotatable_between
            _conf_score = score_pharmacophore_conformance(
                canon, pharm_a_idx, pharm_b_idx, target_dist, tolerance=1.5, n_confs=100)
            _mol_for_rot = Chem.MolFromSmiles(canon)
            _rot_between = count_rotatable_between(_mol_for_rot, pharm_a_idx, pharm_b_idx)
            print(f"  適合率: {_conf_score['conformance_rate']:.1%} "
                  f"(中央値={_conf_score['median_dist']:.1f}Å, "
                  f"回転={_rot_between})")
        except ImportError:
            _conf_score = {"conformance_rate": 0.0, "median_dist": 0.0, "dist_std": 0.0}
            _rot_between = -1

        # ── ⑤ Drug-likeness 計算 ──
        mol_no_h = Chem.RemoveHs(mol_h)
        dl = calc_drug_likeness(mol_no_h)
        print(f"  MW={dl['MW']}, LogP={dl['LogP']}, HBD={dl['HBD']}, "
              f"HBA={dl['HBA']}, PSA={dl['PSA']}, Rot={dl['RotBonds']}")
        print(f"  DrugLike: {dl['DrugLike']} (Ro5={dl['Ro5']}, Veber={dl['Veber']})")

        # ── ⑤b 合成容易性・構造アラート ──
        try:
            from synthesizability import calc_sa_score, calc_qed, check_pains, check_brenk
            _sa = calc_sa_score(mol_no_h)
            _qed = calc_qed(mol_no_h)
            _pains_ok, _pains_alerts = check_pains(mol_no_h)
            _brenk_ok, _brenk_alerts = check_brenk(mol_no_h)
            import math as _math
            _sa_s = f"{_sa:.2f}" if not _math.isnan(_sa) else "N/A"
            _qed_s = f"{_qed:.3f}" if not _math.isnan(_qed) else "N/A"
            print(f"  SA={_sa_s}, QED={_qed_s}, "
                  f"PAINS={'OK' if _pains_ok else 'NG'}, "
                  f"BRENK={'OK' if _brenk_ok else 'NG'}")
        except ImportError:
            _sa = float("nan")
            _qed = float("nan")
            _pains_ok, _pains_alerts = True, []
            _brenk_ok, _brenk_alerts = True, []

        # ── ⑥ SDF 保存 ──
        sdf_path = str(OUT_DIR / f"input_{name}.sdf")
        import math as _math
        props = {
            "AnchorA": anchor_a["name"],
            "Linker":  linker["name"],
            "AnchorB": anchor_b["name"],
            "Res1":    f"{resname1}{resnum1}",
            "Res2":    f"{resname2}{resnum2}",
            "TargetDist": f"{target_dist:.2f}",
            "ActualDist": f"{actual_dist:.2f}",
            **dl,
            "SA_Score":  f"{_sa:.2f}" if not _math.isnan(_sa) else "",
            "QED":       f"{_qed:.4f}" if not _math.isnan(_qed) else "",
            "PAINS_OK":  str(_pains_ok),
            "BRENK_OK":  str(_brenk_ok),
        }
        save_sdf(mol_h, best_cid, name, props, sdf_path)

        # ── ⑦ smina ドッキング ──
        dock_score = None
        docked_sdf = None
        if dock and os.path.isfile(RECEPTOR_PATH) and os.path.isfile(PEPTIDE_REF):
            print(f"  → ドッキング中...")
            dock_score, docked_sdf = run_docking(
                sdf_path, RECEPTOR_PATH, PEPTIDE_REF,
                str(OUT_DIR), name, exhaustiveness,
            )
            if dock_score is not None:
                print(f"  ドッキングスコア: {dock_score:.3f} kcal/mol")
            else:
                print(f"  ドッキング失敗またはスコア取得不可")

        results.append({
            "name":        name,
            "smiles":      canon,
            "res1":        f"{resname1}{resnum1}",
            "res2":        f"{resname2}{resnum2}",
            "anchor_a":   anchor_a["name"],
            "linker":      linker["name"],
            "anchor_b":   anchor_b["name"],
            "target_dist": round(target_dist, 3),
            "actual_dist": round(actual_dist, 3),
            "dock_score":  dock_score,
            "sdf_input":   sdf_path,
            "sdf_docked":  docked_sdf,
            "conformance_rate":  _conf_score["conformance_rate"],
            "median_dist":       _conf_score["median_dist"],
            "dist_std":          _conf_score["dist_std"],
            "rotatable_between": _rot_between,
            **dl,
            "SA_Score":      round(_sa, 2) if not _math.isnan(_sa) else None,
            "QED":           round(_qed, 4) if not _math.isnan(_qed) else None,
            "PAINS_OK":      _pains_ok,
            "BRENK_OK":      _brenk_ok,
            "PAINS_Alerts":  "; ".join(_pains_alerts),
            "BRENK_Alerts":  "; ".join(_brenk_alerts),
        })

    # ── ⑧ 結果をソート・保存 ──
    results.sort(key=lambda r: (
        r["dock_score"] if r["dock_score"] is not None else 0.0,
        -int(r.get("DrugLike", False)),
    ))

    csv_path = OUT_DIR / f"bridge_{res1_str}_{res2_str}_results.csv"
    if results:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\n結果CSV: {csv_path}")

    # ── ⑨ 可視化 ──
    if results:
        _plot_results(results, res1_str, res2_str, target_dist)

    print(f"\n{'='*60}")
    print(f"生成完了: {len(results)} 分子")
    if results and results[0]["dock_score"] is not None:
        best = results[0]
        print(f"最良スコア: {best['dock_score']:.3f} kcal/mol  ({best['name']})")
    print(f"{'='*60}\n")

    return results


# ══════════════════════════════════════════
# 可視化
# ══════════════════════════════════════════

def _plot_results(results: list[dict], res1: str, res2: str, target_dist: float):
    """ドッキングスコア棒グラフと距離誤差散布図を出力"""
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle(f"ファーマコフォアブリッジ: {res1} – {res2}  (目標 {target_dist:.1f} Å)",
                 fontsize=13)

    # ── 左: ドッキングスコア ──
    ax = axes[0]
    scored = [r for r in results if r["dock_score"] is not None]
    if scored:
        names  = [r["name"].split("_", 3)[-1][:30] for r in scored[:15]]
        scores = [r["dock_score"] for r in scored[:15]]
        colors = ["#e74c3c" if r["DrugLike"] else "#95a5a6" for r in scored[:15]]
        bars = ax.barh(names, scores, color=colors)
        ax.set_xlabel("ドッキングスコア (kcal/mol)")
        ax.set_title("ドッキングスコア上位")
        ax.invert_yaxis()
        # 凡例
        patches = [
            mpatches.Patch(color="#e74c3c", label="DrugLike ✓"),
            mpatches.Patch(color="#95a5a6", label="DrugLike ✗"),
        ]
        ax.legend(handles=patches, loc="lower right", fontsize=8)
    else:
        ax.text(0.5, 0.5, "ドッキング未実行", ha="center", va="center",
                transform=ax.transAxes, fontsize=12)
        ax.set_title("ドッキングスコア")

    # ── 右: 実測距離 vs 目標距離 ──
    ax2 = axes[1]
    actual = [r["actual_dist"] for r in results]
    mw_arr  = [r["MW"]  for r in results]
    dl_arr  = [r["DrugLike"] for r in results]
    colors2 = ["#e74c3c" if d else "#95a5a6" for d in dl_arr]

    sc = ax2.scatter(actual, mw_arr, c=colors2, alpha=0.7, s=60, edgecolors="k", linewidths=0.3)
    ax2.axvline(target_dist, color="navy", linestyle="--", linewidth=1.5,
                label=f"目標距離 {target_dist:.1f} Å")
    ax2.axvspan(target_dist - 2.5, target_dist + 2.5, alpha=0.1, color="navy")
    ax2.set_xlabel("実測ファーマコフォア間距離 (Å)")
    ax2.set_ylabel("分子量 (Da)")
    ax2.set_title("距離充足 vs 分子量")
    ax2.legend(fontsize=8)

    plt.tight_layout()
    out_png = OUT_DIR / f"bridge_{res1}_{res2}_analysis.png"
    plt.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"グラフ保存: {out_png}")


def plot_pharmacophore_geometry(
    res1_str: str, res2_str: str,
    pdb_path: str = PDB_PATH, chain_id: str = PROTEIN_CHAIN,
):
    """ファーマコフォア点とCβ位置の3D俯瞰図を出力 (参考用)"""
    try:
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
    except ImportError:
        return

    resname1, resnum1 = parse_residue_id(res1_str)
    resname2, resnum2 = parse_residue_id(res2_str)
    coords1, _ = get_cbeta_coords(pdb_path, chain_id, resnum1)
    coords2, _ = get_cbeta_coords(pdb_path, chain_id, resnum2)
    dist = np.linalg.norm(coords1 - coords2)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(6, 5))
    ax  = fig.add_subplot(111, projection="3d")
    ax.scatter(*coords1, c="red",  s=200, label=f"{res1_str} Cβ", zorder=5)
    ax.scatter(*coords2, c="blue", s=200, label=f"{res2_str} Cβ", zorder=5)
    ax.plot([coords1[0], coords2[0]],
            [coords1[1], coords2[1]],
            [coords1[2], coords2[2]], "k--", linewidth=2)
    mid = (coords1 + coords2) / 2
    ax.text(mid[0], mid[1], mid[2], f"{dist:.1f} Å", fontsize=10, ha="center")
    ax.set_title(f"Cβ–Cβ 距離: {res1_str} — {res2_str}")
    ax.legend()
    plt.tight_layout()
    out = OUT_DIR / f"geometry_{res1_str}_{res2_str}.png"
    plt.savefig(out, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"ジオメトリ図: {out}")


# ══════════════════════════════════════════
# エントリポイント
# ══════════════════════════════════════════

def main():
    parser = argparse.ArgumentParser(
        description="ファーマコフォアブリッジ低分子生成",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  python pharmacophore_bridge.py --point1 LYS48 --point2 LEU50
  python pharmacophore_bridge.py --point1 LYS46 --point2 TYR49 --no-dock
  python pharmacophore_bridge.py --point1 ARG28 --point2 TYR49 --n-confs 100
  python pharmacophore_bridge.py --list-pairs
        """,
    )
    parser.add_argument("--point1",  type=str, help="残基1 (例: LYS48)")
    parser.add_argument("--point2",  type=str, help="残基2 (例: LEU50)")
    parser.add_argument("--chain",   type=str, default=PROTEIN_CHAIN,
                        help="タンパク質チェーンID (デフォルト: A)")
    parser.add_argument("--no-dock", action="store_true",
                        help="ドッキングをスキップ")
    parser.add_argument("--n-confs", type=int, default=50,
                        help="コンフォマー生成数 (デフォルト: 50)")
    parser.add_argument("--tolerance", type=float, default=2.5,
                        help="距離制約の許容範囲 Å (デフォルト: 2.5)")
    parser.add_argument("--exhaustiveness", type=int, default=8,
                        help="smina 探索徹底度 (デフォルト: 8)")
    parser.add_argument("--list-pairs", action="store_true",
                        help="既知ポケットペアを表示して終了")
    parser.add_argument("--all-pairs", action="store_true",
                        help="全既知ペアを順次処理")
    args = parser.parse_args()

    if args.list_pairs:
        print("\n既知ポケットペア一覧:")
        print(f"{'残基1':10} {'残基2':10} {'Cβ–Cβ距離':>12} {'タイプ':15}")
        print("-" * 52)
        for p in KNOWN_PAIRS:
            print(f"{p['res1']:10} {p['res2']:10} {p['dist']:>10.1f} Å  {p['type']}")
        return

    if args.all_pairs:
        for pair in KNOWN_PAIRS:
            generate_bridges(
                pair["res1"], pair["res2"],
                chain_id=args.chain,
                dock=not args.no_dock,
                n_confs=args.n_confs,
                tolerance=args.tolerance,
                exhaustiveness=args.exhaustiveness,
            )
        return

    if not args.point1 or not args.point2:
        parser.print_help()
        print("\n⚠ --point1 と --point2 を指定してください。")
        print("  既知ペアを確認: python pharmacophore_bridge.py --list-pairs")
        sys.exit(1)

    plot_pharmacophore_geometry(args.point1, args.point2,
                                chain_id=args.chain)
    results = generate_bridges(
        args.point1, args.point2,
        chain_id=args.chain,
        dock=not args.no_dock,
        n_confs=args.n_confs,
        tolerance=args.tolerance,
        exhaustiveness=args.exhaustiveness,
    )

    if results:
        print("\n上位5候補:")
        print(f"{'順位':4} {'スコア':>10} {'DrugLike':>9} {'実測dist':>10} {'名前'}")
        print("-" * 65)
        for i, r in enumerate(results[:5], 1):
            sc = f"{r['dock_score']:.3f}" if r['dock_score'] is not None else "  N/A "
            dl = "✓" if r["DrugLike"] else "✗"
            print(f"{i:4}  {sc:>10}  {dl:>9}  {r['actual_dist']:>8.2f} Å  {r['name']}")


if __name__ == "__main__":
    main()
