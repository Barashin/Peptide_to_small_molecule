"""
rigid_scaffold_design.py
========================
剛直スキャフォールドを用いたファーマコフォアブリッジ設計 + コンフォメーション適合性スコアリング

アプローチ A: 既知の剛直骨格 (ベンゼン, ナフタレン, インドール, ノルボルナン, ピペラジン 等)
            を置換基間距離でフィルタし, アンカーA + スキャフォールド + アンカーB を組立て
アプローチ B: 生成された分子の低エネルギーコンフォマーを多数生成し,
            ファーマコフォア距離の維持率 (conformance_rate) と
            ファーマコフォア間の回転可能結合数 (rotatable_between) を評価

使い方:
  python rigid_scaffold_design.py --point1 LYS48 --point2 LEU50
  python rigid_scaffold_design.py --point1 LYS46 --point2 TYR49 --no-dock
  python rigid_scaffold_design.py --list-scaffolds
"""

from __future__ import annotations

import csv
import math
import os
import re
import statistics
import sys
from pathlib import Path

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdDistGeom, rdMolDescriptors
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ──────────────────────────────────────────
# pharmacophore_bridge から共通関数をインポート
# ──────────────────────────────────────────
from pharmacophore_bridge import (
    ANCHOR_A,
    ANCHOR_B,
    PDB_PATH,
    PROTEIN_CHAIN,
    RECEPTOR_PATH,
    PEPTIDE_REF,
    parse_residue_id,
    get_cbeta_coords,
    embed_with_constraint,
    select_best_conformer,
    save_sdf,
    run_docking,
)
from utils.residue_defs import RES_TYPE
from utils.drug_likeness import calculate_drug_likeness
import linker_library as _ll

# ──────────────────────────────────────────
# パス設定
# ──────────────────────────────────────────
BASE_DIR = Path(__file__).parent
OUT_DIR  = BASE_DIR / "results" / "rigid_scaffold"

# ══════════════════════════════════════════
# 剛直スキャフォールドライブラリ
# ══════════════════════════════════════════
#
# 各エントリ:
#   name        : スキャフォールド名
#   smiles_r1r2 : [*:1] と [*:2] で置換位置を示す SMILES
#   dist_r1r2   : 置換基間距離 (Å, 3D 最適化で事前計算)
#   n_rotatable : コア内の回転可能結合数
#   ring_type   : aromatic / saturated / mixed
#   desc        : 説明
#
# 距離は 隣接 C–C 結合 ~1.5 Å (飽和) / ~1.4 Å (芳香族) を基準に
# 各置換パターンの 3D 距離をRDKit ETKDG + MMFF で測定した値。

RIGID_SCAFFOLDS: list[dict] = [
    # ── 単環芳香族 ──
    {
        "name": "benzene_ortho",
        "smiles_r1r2": "[*:1]c1ccccc1[*:2]",
        "dist_r1r2": 2.8,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ベンゼン ortho (1,2-二置換)",
    },
    {
        "name": "benzene_meta",
        "smiles_r1r2": "[*:1]c1cccc([*:2])c1",
        "dist_r1r2": 4.3,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ベンゼン meta (1,3-二置換)",
    },
    {
        "name": "benzene_para",
        "smiles_r1r2": "[*:1]c1ccc([*:2])cc1",
        "dist_r1r2": 5.7,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ベンゼン para (1,4-二置換)",
    },
    {
        "name": "pyridine_2_4",
        "smiles_r1r2": "[*:1]c1cc([*:2])ccn1",
        "dist_r1r2": 4.3,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ピリジン 2,4-二置換",
    },
    {
        "name": "pyridine_2_5",
        "smiles_r1r2": "[*:1]c1ccc([*:2])cn1",
        "dist_r1r2": 5.5,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ピリジン 2,5-二置換",
    },
    {
        "name": "pyrimidine_2_5",
        "smiles_r1r2": "[*:1]c1ncc([*:2])cn1",
        "dist_r1r2": 5.5,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ピリミジン 2,5-二置換",
    },

    # ── 縮環芳香族 ──
    {
        "name": "naphthalene_1_4",
        "smiles_r1r2": "[*:1]c1ccc([*:2])c2ccccc12",
        "dist_r1r2": 5.0,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ナフタレン 1,4-二置換",
    },
    {
        "name": "naphthalene_2_6",
        "smiles_r1r2": "[*:1]c1ccc2cc([*:2])ccc2c1",
        "dist_r1r2": 8.6,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ナフタレン 2,6-二置換",
    },
    {
        "name": "naphthalene_1_5",
        "smiles_r1r2": "[*:1]c1cccc2c([*:2])cccc12",
        "dist_r1r2": 6.9,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ナフタレン 1,5-二置換",
    },
    {
        "name": "indole_3_5",
        "smiles_r1r2": "[*:1]c1c[nH]c2ccc([*:2])cc12",
        "dist_r1r2": 5.3,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "インドール 3,5-二置換",
    },
    {
        "name": "indole_3_6",
        "smiles_r1r2": "[*:1]c1c[nH]c2cc([*:2])ccc12",
        "dist_r1r2": 6.3,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "インドール 3,6-二置換",
    },
    {
        "name": "benzimidazole_2_5",
        "smiles_r1r2": "[*:1]c1nc2ccc([*:2])cc2[nH]1",
        "dist_r1r2": 5.5,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ベンズイミダゾール 2,5-二置換",
    },
    {
        "name": "quinoline_2_6",
        "smiles_r1r2": "[*:1]c1ccc2cc([*:2])ccc2n1",
        "dist_r1r2": 8.5,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "キノリン 2,6-二置換",
    },
    {
        "name": "quinoline_2_7",
        "smiles_r1r2": "[*:1]c1ccc2c([*:2])cccc2n1",
        "dist_r1r2": 7.0,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "キノリン 2,7-二置換",
    },

    # ── 飽和二環 (高い剛直性) ──
    {
        "name": "norbornane_2_5",
        "smiles_r1r2": "[*:1]C1CC2CC1CC2[*:2]",
        "dist_r1r2": 3.9,
        "n_rotatable": 0,
        "ring_type": "saturated",
        "desc": "ノルボルナン 2,5-二置換",
    },
    {
        "name": "BCO_1_4",
        "smiles_r1r2": "[*:1]C1CCC([*:2])CC1",
        "dist_r1r2": 3.9,
        "n_rotatable": 0,
        "ring_type": "saturated",
        "desc": "ビシクロ[2.2.2]オクタン 1,4-二置換",
    },
    {
        "name": "adamantane_1_3",
        "smiles_r1r2": "[*:1]C12CC3CC(CC([*:2])C3)C1C2",
        "dist_r1r2": 3.7,
        "n_rotatable": 0,
        "ring_type": "saturated",
        "desc": "アダマンタン 1,3-二置換",
    },
    # ── 含窒素ヘテロ環 (飽和) ──
    {
        "name": "piperazine_1_4",
        "smiles_r1r2": "[*:1]N1CCN([*:2])CC1",
        "dist_r1r2": 3.9,
        "n_rotatable": 0,
        "ring_type": "saturated",
        "desc": "ピペラジン 1,4-二置換",
    },
    {
        "name": "piperidine_1_4",
        "smiles_r1r2": "[*:1]N1CCC([*:2])CC1",
        "dist_r1r2": 3.9,
        "n_rotatable": 0,
        "ring_type": "saturated",
        "desc": "ピペリジン 1,4-二置換",
    },
    {
        "name": "morpholine_2_5",
        "smiles_r1r2": "[*:1]C1CN([*:2])CCO1",
        "dist_r1r2": 3.9,
        "n_rotatable": 0,
        "ring_type": "saturated",
        "desc": "モルホリン 2,5-二置換",
    },

    # ── sp2 リンカー含有 (やや長距離) ──
    {
        "name": "biphenyl_4_4",
        "smiles_r1r2": "[*:1]c1ccc(-c2ccc([*:2])cc2)cc1",
        "dist_r1r2": 10.2,
        "n_rotatable": 1,
        "ring_type": "aromatic",
        "desc": "ビフェニル 4,4'-二置換",
    },
    {
        "name": "trans_stilbene",
        "smiles_r1r2": "[*:1]c1ccc(/C=C/c2ccc([*:2])cc2)cc1",
        "dist_r1r2": 12.6,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "トランススチルベン 4,4'-二置換",
    },

    # ── 縮環ヘテロ環 ──
    {
        "name": "benzothiazole_2_6",
        "smiles_r1r2": "[*:1]c1nc2ccc([*:2])cc2s1",
        "dist_r1r2": 6.6,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ベンゾチアゾール 2,6-二置換",
    },
    {
        "name": "benzoxazole_2_5",
        "smiles_r1r2": "[*:1]c1nc2ccc([*:2])cc2o1",
        "dist_r1r2": 5.5,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "ベンゾオキサゾール 2,5-二置換",
    },
    {
        "name": "isoquinoline_1_6",
        "smiles_r1r2": "[*:1]c1nccc2cc([*:2])ccc12",
        "dist_r1r2": 7.5,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "イソキノリン 1,6-二置換",
    },
    {
        "name": "imidazopyridine_2_6",
        "smiles_r1r2": "[*:1]c1nc2ccc([*:2])cn2c1",
        "dist_r1r2": 6.0,
        "n_rotatable": 0,
        "ring_type": "aromatic",
        "desc": "イミダゾ[1,2-a]ピリジン 2,6-二置換",
    },
]

# スキャフォールド距離の検証 & キャッシュ (初回ロード時に実行)
_SCAFFOLDS_VALIDATED = False


def _validate_scaffold_distances():
    """RIGID_SCAFFOLDS の dist_r1r2 を RDKit 3D で検証・補正する (初回のみ)"""
    global _SCAFFOLDS_VALIDATED
    if _SCAFFOLDS_VALIDATED:
        return

    from rdkit import RDLogger
    logger = RDLogger.logger()
    logger.setLevel(RDLogger.ERROR)  # UFF 警告を抑制

    for sc in RIGID_SCAFFOLDS:
        smi = sc["smiles_r1r2"]
        # ダミー原子 [*:1], [*:2] をメチル基に置換して 3D 埋め込み
        test_smi = smi.replace("[*:1]", "[CH3:1]").replace("[*:2]", "[CH3:2]")
        mol = Chem.MolFromSmiles(test_smi)
        if mol is None:
            continue

        # アトムマップ 1, 2 のメチル C のインデックスを取得
        idx1 = idx2 = None
        for atom in mol.GetAtoms():
            mp = atom.GetAtomMapNum()
            if mp == 1:
                idx1 = atom.GetIdx()
                atom.SetAtomMapNum(0)
            elif mp == 2:
                idx2 = atom.GetIdx()
                atom.SetAtomMapNum(0)
        if idx1 is None or idx2 is None:
            continue

        mol_h = Chem.AddHs(mol)
        try:
            params = AllChem.ETKDGv3()
            params.randomSeed = 42
            cid = AllChem.EmbedMolecule(mol_h, params)
            if cid < 0:
                continue
            AllChem.MMFFOptimizeMolecule(mol_h, confId=cid, maxIters=500)
        except Exception:
            continue

        conf = mol_h.GetConformer(cid)
        pos1 = np.array(conf.GetAtomPosition(idx1))
        pos2 = np.array(conf.GetAtomPosition(idx2))
        measured = float(np.linalg.norm(pos1 - pos2))
        sc["dist_r1r2"] = round(measured, 1)

    logger.setLevel(RDLogger.WARNING)
    _SCAFFOLDS_VALIDATED = True


# ══════════════════════════════════════════
# スキャフォールド選択
# ══════════════════════════════════════════

def select_scaffolds_by_dist(
    target_dist: float,
    anchor_a_reach: float,
    anchor_b_reach: float,
    tolerance: float = 1.5,
    max_n: int = 8,
) -> list[dict]:
    """
    アンカー到達距離を差し引いた残距離にマッチする剛直スキャフォールドを選択する。

    Args:
        target_dist     : Cβ–Cβ 間距離 (Å)
        anchor_a_reach  : アンカー A の到達距離 (Å)
        anchor_b_reach  : アンカー B の到達距離 (Å)
        tolerance       : 許容誤差 (Å)
        max_n           : 最大返却数

    Returns:
        マッチしたスキャフォールドの list[dict]
    """
    _validate_scaffold_distances()

    scaffold_needed = target_dist - anchor_a_reach - anchor_b_reach

    hits = []
    for sc in RIGID_SCAFFOLDS:
        dist_err = abs(sc["dist_r1r2"] - scaffold_needed)
        if dist_err <= tolerance:
            hits.append((dist_err, sc["n_rotatable"], sc))

    # 距離誤差の昇順 → 回転可能結合数の昇順でソート
    hits.sort(key=lambda x: (x[0], x[1]))

    result = [h[2] for h in hits[:max_n]]

    # ヒットなしの場合: 最近傍のスキャフォールドを1つ返す
    if not result:
        all_sorted = sorted(
            RIGID_SCAFFOLDS,
            key=lambda sc: abs(sc["dist_r1r2"] - scaffold_needed),
        )
        if all_sorted:
            result = [all_sorted[0]]

    return result


# ══════════════════════════════════════════
# SMILES 組み立て
# ══════════════════════════════════════════

def _find_connect_atom(mol: Chem.Mol, pharm_offset: int) -> int | None:
    """
    アンカー分子内でスキャフォールドと接続すべき原子を見つける。
    pharm_offset (ファーマコフォア原子) から最も遠い、implicit H を持つ原子を返す。
    (implicit H は新しい結合で置換できる)
    """
    n = mol.GetNumAtoms()
    if n == 1:
        return 0  # 単原子の場合 (例: N, O)

    candidates = []
    for i in range(n):
        atom = mol.GetAtomWithIdx(i)
        if atom.GetAtomicNum() == 0:
            continue
        # implicit H を持つ原子は新しい結合を受け入れ可能
        if atom.GetTotalNumHs() > 0:
            if i == pharm_offset:
                dist = 0  # 同一原子 (単原子アンカーでは唯一の選択肢)
            else:
                path = Chem.rdmolops.GetShortestPath(mol, pharm_offset, i)
                dist = len(path) - 1
            candidates.append((dist, i))

    if not candidates:
        return None
    # ファーマコフォアから最も遠い原子を選択
    candidates.sort(key=lambda x: -x[0])
    return candidates[0][1]


def build_rigid_bridge_smiles(
    anchor_a: dict,
    scaffold: dict,
    anchor_b: dict,
) -> tuple[str, int, int] | None:
    """
    アンカー A + 剛直スキャフォールド + アンカー B の SMILES を組み立て、
    ファーマコフォア原子インデックスを返す。

    RDKit の CombineMols + RWMol.AddBond でスキャフォールドの [*:1], [*:2]
    ダミー原子を除去してアンカーを接続する。
    アトムマップ番号 91/92 でファーマコフォア原子を追跡する。

    Returns:
        (canonical_smiles, pharm_a_idx, pharm_b_idx) または None
    """
    # スキャフォールドの解析
    mol_scaffold = Chem.MolFromSmiles(scaffold["smiles_r1r2"])
    if mol_scaffold is None:
        return None

    # ダミー原子 [*:1], [*:2] のインデックスとその隣接原子を取得
    dummy1_idx = dummy2_idx = None
    attach1_idx = attach2_idx = None
    for atom in mol_scaffold.GetAtoms():
        mp = atom.GetAtomMapNum()
        if mp == 1:
            dummy1_idx = atom.GetIdx()
            neighbors = atom.GetNeighbors()
            if neighbors:
                attach1_idx = neighbors[0].GetIdx()
        elif mp == 2:
            dummy2_idx = atom.GetIdx()
            neighbors = atom.GetNeighbors()
            if neighbors:
                attach2_idx = neighbors[0].GetIdx()

    if None in (dummy1_idx, dummy2_idx, attach1_idx, attach2_idx):
        return None

    # アンカー Mol を生成
    mol_a = Chem.MolFromSmiles(anchor_a["smiles"])
    mol_b = Chem.MolFromSmiles(anchor_b["smiles"])
    if mol_a is None or mol_b is None:
        return None

    pharm_offset_a = anchor_a["pharm_offset"]
    pharm_offset_b = anchor_b["pharm_offset"]
    n_a = mol_a.GetNumAtoms()
    n_b = mol_b.GetNumAtoms()

    # アンカーの接続原子を決定 (利用可能な原子価を持ち、ファーマコフォアから最遠)
    connect_a_idx = _find_connect_atom(mol_a, pharm_offset_a)
    connect_b_idx = _find_connect_atom(mol_b, pharm_offset_b)
    if connect_a_idx is None or connect_b_idx is None:
        return None

    # ファーマコフォア原子にアトムマップを設定
    mol_a = Chem.RWMol(mol_a)
    if pharm_offset_a < n_a:
        mol_a.GetAtomWithIdx(pharm_offset_a).SetAtomMapNum(91)
    mol_b = Chem.RWMol(mol_b)
    if pharm_offset_b < n_b:
        mol_b.GetAtomWithIdx(pharm_offset_b).SetAtomMapNum(92)

    # CombineMols で結合: scaffold + anchor_a + anchor_b
    combined = Chem.CombineMols(mol_scaffold, mol_a)
    combined = Chem.CombineMols(combined, mol_b)
    combined_rw = Chem.RWMol(combined)

    # インデックスのオフセット計算
    n_scaffold = mol_scaffold.GetNumAtoms()
    offset_a = n_scaffold
    offset_b = n_scaffold + n_a

    # スキャフォールド接続点とアンカー接続点に結合を追加
    combined_rw.AddBond(attach1_idx, offset_a + connect_a_idx,
                        Chem.BondType.SINGLE)
    combined_rw.AddBond(attach2_idx, offset_b + connect_b_idx,
                        Chem.BondType.SINGLE)

    # ダミー原子を除去 (インデックスが大きい方から除去)
    to_remove = sorted([dummy1_idx, dummy2_idx], reverse=True)
    for idx in to_remove:
        combined_rw.RemoveAtom(idx)

    # サニタイズ
    try:
        Chem.SanitizeMol(combined_rw)
    except Exception:
        return None

    # アトムマップ付き SMILES を取得 → 再パースして正規インデックスを得る
    mapped_smi = Chem.MolToSmiles(combined_rw)
    mol_mapped = Chem.MolFromSmiles(mapped_smi)
    if mol_mapped is None:
        return None

    # マップ付きSMILES内でのインデックスを取得
    final_a = final_b = None
    for atom in mol_mapped.GetAtoms():
        mp = atom.GetAtomMapNum()
        if mp == 91:
            final_a = atom.GetIdx()
        elif mp == 92:
            final_b = atom.GetIdx()

    if final_a is None or final_b is None:
        return None

    # マップを除去してクリーンSMILESを取得
    # 重要: MolToSmiles でマップなしSMILESを取得すると
    #        原子順序が変わるため、マップ付きSMILESのまま返却し
    #        呼び出し元でマップを除去する
    return mapped_smi, final_a, final_b


# ══════════════════════════════════════════
# コンフォメーション適合性スコアリング (アプローチ B)
# ══════════════════════════════════════════

def score_pharmacophore_conformance(
    smiles: str,
    pharm_a_idx: int,
    pharm_b_idx: int,
    target_dist: float,
    tolerance: float = 1.5,
    n_confs: int = 200,
    energy_window: float = 5.0,
) -> dict:
    """
    制約なしコンフォマー探索でファーマコフォア距離の維持率を評価する。

    Args:
        smiles         : 分子の SMILES
        pharm_a_idx    : ファーマコフォア原子 A のインデックス
        pharm_b_idx    : ファーマコフォア原子 B のインデックス
        target_dist    : 目標ファーマコフォア距離 (Å)
        tolerance      : 許容誤差 (Å)
        n_confs        : 生成するコンフォマー数
        energy_window  : エネルギーウィンドウ (kcal/mol)

    Returns:
        {
            "conformance_rate": float,  # 目標距離内のコンフォマー割合 (0–1)
            "median_dist": float,       # 中央値距離 (Å)
            "dist_std": float,          # 距離の標準偏差 (Å)
            "n_confs_valid": int,       # 有効コンフォマー数
        }
    """
    default = {
        "conformance_rate": 0.0,
        "median_dist": 0.0,
        "dist_std": 0.0,
        "n_confs_valid": 0,
    }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return default

    mol_h = Chem.AddHs(mol)

    # 制約なしコンフォマー生成
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    params.numThreads = 0  # 全コアを使用
    conf_ids = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params)
    if not conf_ids:
        return default

    # MMFF 最適化
    opt_results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, maxIters=500, numThreads=0)
    if not opt_results:
        return default

    # エネルギーフィルタ
    energies = []
    for cid, (converged, energy) in zip(conf_ids, opt_results):
        if converged == -1:  # 最適化失敗
            continue
        energies.append((cid, energy))

    if not energies:
        return default

    min_energy = min(e for _, e in energies)
    valid_confs = [(cid, e) for cid, e in energies
                   if e <= min_energy + energy_window]

    if not valid_confs:
        return default

    # 各コンフォマーのファーマコフォア距離を測定
    distances = []
    for cid, _ in valid_confs:
        conf = mol_h.GetConformer(cid)
        pos_a = np.array(conf.GetAtomPosition(pharm_a_idx))
        pos_b = np.array(conf.GetAtomPosition(pharm_b_idx))
        dist = float(np.linalg.norm(pos_a - pos_b))
        distances.append(dist)

    # 適合率計算
    n_conforming = sum(
        1 for d in distances
        if target_dist - tolerance <= d <= target_dist + tolerance
    )
    conformance_rate = n_conforming / len(distances)

    return {
        "conformance_rate": round(conformance_rate, 4),
        "median_dist": round(statistics.median(distances), 2),
        "dist_std": round(statistics.stdev(distances), 2) if len(distances) > 1 else 0.0,
        "n_confs_valid": len(valid_confs),
    }


def count_rotatable_between(mol: Chem.Mol, idx_a: int, idx_b: int) -> int:
    """
    ファーマコフォア原子 A–B 間の最短パス上の回転可能結合数を数える。

    環内結合は RDKit の Bond.GetIsRotatable() で自動的に除外される。

    Args:
        mol   : RDKit Mol オブジェクト (水素なし)
        idx_a : ファーマコフォア原子 A のインデックス
        idx_b : ファーマコフォア原子 B のインデックス

    Returns:
        回転可能結合数 (int)
    """
    if mol is None:
        return -1

    try:
        path = Chem.rdmolops.GetShortestPath(mol, idx_a, idx_b)
    except Exception:
        return -1

    if len(path) < 2:
        return 0

    n_rot = 0
    for i in range(len(path) - 1):
        bond = mol.GetBondBetweenAtoms(path[i], path[i + 1])
        if bond is None:
            continue
        # 環内結合は回転可能ではない
        if bond.IsInRing():
            continue
        # 単結合かつ回転可能
        if bond.GetBondTypeAsDouble() == 1.0 and not bond.IsInRing():
            n_rot += 1

    return n_rot


# ══════════════════════════════════════════
# 既存ブリッジ結果への適合性スコア追加
# ══════════════════════════════════════════

def add_conformance_to_existing(
    bridge_results: list[dict],
    target_dist: float,
    tolerance: float = 1.5,
) -> list[dict]:
    """
    既存の pharmacophore_bridge.generate_bridges() の結果に
    conformance_rate と rotatable_between を追加する。

    Args:
        bridge_results : generate_bridges() の返却リスト
        target_dist    : 目標ファーマコフォア距離 (Å)
        tolerance      : 適合判定の許容範囲 (Å)

    Returns:
        conformance_rate, rotatable_between が追加された bridge_results
    """
    for rec in bridge_results:
        smiles = rec.get("smiles", "")
        if not smiles:
            rec["conformance_rate"] = 0.0
            rec["median_dist"] = 0.0
            rec["dist_std"] = 0.0
            rec["rotatable_between"] = -1
            continue

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            rec["conformance_rate"] = 0.0
            rec["median_dist"] = 0.0
            rec["dist_std"] = 0.0
            rec["rotatable_between"] = -1
            continue

        # ファーマコフォアインデックスの推定
        # (bridge_results の anchor 情報から再構築)
        built = _rebuild_pharm_indices(rec, mol)
        if built is None:
            # フォールバック: 最も離れた重原子ペアを使用
            pharm_a, pharm_b = _find_farthest_atoms(mol)
        else:
            pharm_a, pharm_b = built

        # コンフォメーション適合性スコア
        conf_score = score_pharmacophore_conformance(
            smiles, pharm_a, pharm_b, target_dist, tolerance,
            n_confs=100,  # 既存結果は軽量スコアリング
        )
        rec["conformance_rate"] = conf_score["conformance_rate"]
        rec["median_dist"] = conf_score["median_dist"]
        rec["dist_std"] = conf_score["dist_std"]

        # 回転可能結合数
        rec["rotatable_between"] = count_rotatable_between(mol, pharm_a, pharm_b)

    return bridge_results


def _rebuild_pharm_indices(rec: dict, mol: Chem.Mol) -> tuple[int, int] | None:
    """結果dictからファーマコフォアインデックスを再構築する (ベストエフォート)"""
    # anchor_a/anchor_b の名前から ANCHOR_A/B の pharm_offset を取得
    anchor_a_name = rec.get("anchor_a", "")
    anchor_b_name = rec.get("anchor_b", "")
    if not anchor_a_name or not anchor_b_name:
        return None

    # pharm_offset を取得
    offset_a = None
    for type_anchors in ANCHOR_A.values():
        for a in type_anchors:
            if a["name"] == anchor_a_name:
                offset_a = a["pharm_offset"]
                break
    offset_b = None
    for type_anchors in ANCHOR_B.values():
        for a in type_anchors:
            if a["name"] == anchor_b_name:
                offset_b = a["pharm_offset"]
                break

    if offset_a is None or offset_b is None:
        return None

    n_atoms = mol.GetNumHeavyAtoms()
    # 簡易推定: A は先頭付近, B は末尾付近
    pharm_a = min(offset_a, n_atoms - 1)
    pharm_b = min(n_atoms - 1 - offset_b, n_atoms - 1)
    if pharm_a == pharm_b:
        pharm_b = n_atoms - 1
    return pharm_a, pharm_b


def _find_farthest_atoms(mol: Chem.Mol) -> tuple[int, int]:
    """3D 座標がない場合のフォールバック: トポロジカルに最も離れた原子ペアを返す"""
    n = mol.GetNumHeavyAtoms()
    if n < 2:
        return 0, max(0, n - 1)
    dm = Chem.rdmolops.GetDistanceMatrix(mol)
    max_dist = 0
    best = (0, n - 1)
    for i in range(n):
        for j in range(i + 1, n):
            if dm[i][j] > max_dist:
                max_dist = dm[i][j]
                best = (i, j)
    return best


# ══════════════════════════════════════════
# メイン生成関数
# ══════════════════════════════════════════

def generate_rigid_bridges(
    res1_str: str,
    res2_str: str,
    pdb_path: str  = PDB_PATH,
    chain_id: str  = PROTEIN_CHAIN,
    dock: bool     = True,
    n_confs: int   = 50,
    tolerance: float = 2.0,
    exhaustiveness: int = 8,
) -> list[dict]:
    """
    2残基間を剛直スキャフォールドで橋渡しする低分子を生成する。

    pharmacophore_bridge.generate_bridges() と同様のフローだが、
    線形リンカーの代わりに RIGID_SCAFFOLDS を使用する。

    Args:
        res1_str : 残基1 (例: "LYS48")
        res2_str : 残基2 (例: "LEU50")
        pdb_path : タンパク質PDBファイルパス
        chain_id : ペプチドチェーンID
        dock     : smina でドッキングするか
        n_confs  : コンフォマー生成数
        tolerance: 距離許容誤差 (Å)
        exhaustiveness: smina exhaustiveness

    Returns:
        results: 各候補分子の情報を含む dict のリスト
    """
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    # ── ① 残基情報・座標取得 ──
    resname1, resnum1 = parse_residue_id(res1_str)
    resname2, resnum2 = parse_residue_id(res2_str)

    print(f"\n{'='*60}")
    print(f"剛直スキャフォールドブリッジ生成")
    print(f"  残基1: {resname1}{resnum1}  残基2: {resname2}{resnum2}")

    coords1, rn1 = get_cbeta_coords(pdb_path, chain_id, resnum1)
    coords2, rn2 = get_cbeta_coords(pdb_path, chain_id, resnum2)
    target_dist = float(np.linalg.norm(coords1 - coords2))
    print(f"  Cβ–Cβ 距離: {target_dist:.2f} Å")

    type1 = RES_TYPE.get(rn1, "HYD")
    type2 = RES_TYPE.get(rn2, "HYD")
    print(f"  残基タイプ: {rn1}({type1}) — {rn2}({type2})")

    # ── ② アンカーとスキャフォールドを選択 ──
    anchors_a = ANCHOR_A.get(type1, ANCHOR_A["HYD"])
    anchors_b = ANCHOR_B.get(type2, ANCHOR_B["HYD"])

    print(f"  アンカーA ({type1}): {[a['name'] for a in anchors_a]}")
    print(f"  アンカーB ({type2}): {[a['name'] for a in anchors_b]}")

    # アンカーペアごとにスキャフォールドを選択
    combos: list[tuple] = []
    scaffold_names_seen: set[str] = set()
    for anchor_a in anchors_a:
        for anchor_b in anchors_b:
            a_reach = _ll.estimate_anchor_reach(anchor_a["smiles"])
            b_reach = _ll.estimate_anchor_reach(anchor_b["smiles"])
            scaffolds = select_scaffolds_by_dist(
                target_dist, a_reach, b_reach, tolerance=tolerance,
            )
            for sc in scaffolds:
                combos.append((anchor_a, sc, anchor_b))
                scaffold_names_seen.add(sc["name"])

    print(f"  スキャフォールド候補: {sorted(scaffold_names_seen)}")
    print(f"  組み合わせ数: {len(combos)}")
    print(f"{'='*60}")

    results  = []
    seen_smi = set()
    mol_idx  = 0

    for anchor_a, scaffold, anchor_b in combos:
        built = build_rigid_bridge_smiles(anchor_a, scaffold, anchor_b)
        if built is None:
            continue
        mapped_smiles, pharm_a_idx, pharm_b_idx = built

        # マップ付き SMILES をパース (インデックスが正しい)
        mol_mapped = Chem.MolFromSmiles(mapped_smiles)
        if mol_mapped is None:
            continue

        # クリーン SMILES (マップなし) を生成して重複排除
        mol_clean = Chem.RWMol(mol_mapped)
        for atom in mol_clean.GetAtoms():
            atom.SetAtomMapNum(0)
        canon = Chem.MolToSmiles(mol_clean)
        if canon in seen_smi:
            continue
        seen_smi.add(canon)

        # Ro5 フィルタ (MWのみ事前確認)
        mw = Descriptors.ExactMolWt(mol_clean)
        if mw > 600:
            continue

        mol_idx += 1
        name = (f"rigid_{res1_str}_{res2_str}"
                f"_{anchor_a['name']}_{scaffold['name']}_{anchor_b['name']}")
        name = name.replace("(", "").replace(")", "").replace(" ", "_")[:60]

        print(f"\n[{mol_idx:03d}] {name}")
        print(f"  SMILES: {canon}")
        print(f"  MW={mw:.1f}, scaffold={scaffold['name']} (dist={scaffold['dist_r1r2']:.1f}Å)")
        print(f"  pharm_A={pharm_a_idx}, pharm_B={pharm_b_idx}")

        # ── ③ 制約付き3D埋め込み (マップ付きSMILESを使用してインデックスの整合性を保証) ──
        mol_h, conf_ids = embed_with_constraint(
            mapped_smiles, pharm_a_idx, pharm_b_idx,
            target_dist, tolerance + 0.5, n_confs, seed=42 + mol_idx,
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

        # ── ⑤ コンフォメーション適合性スコア (アプローチ B) ──
        conf_score = score_pharmacophore_conformance(
            mapped_smiles, pharm_a_idx, pharm_b_idx, target_dist, tolerance=1.5,
        )
        print(f"  適合率: {conf_score['conformance_rate']:.1%} "
              f"(中央値={conf_score['median_dist']:.1f}Å, σ={conf_score['dist_std']:.2f}Å, "
              f"n={conf_score['n_confs_valid']})")

        # ── ⑥ 回転可能結合数 (ファーマコフォア間) ──
        mol_no_h = Chem.RemoveHs(mol_h)
        rot_between = count_rotatable_between(mol_no_h, pharm_a_idx, pharm_b_idx)
        print(f"  ファーマコフォア間回転可能結合: {rot_between}")

        # ── ⑦ Drug-likeness 計算 ──
        dl = calculate_drug_likeness(mol_no_h)
        print(f"  MW={dl['MW']}, LogP={dl['LogP']}, HBD={dl['HBD']}, "
              f"HBA={dl['HBA']}, PSA={dl['PSA']}, Rot={dl['RotBonds']}")
        print(f"  DrugLike: {dl['DrugLike']} (Ro5={dl['Ro5']}, Veber={dl['Veber']})")

        # ── ⑧ 合成容易性・構造アラート ──
        try:
            from synthesizability import calc_sa_score, calc_qed, check_pains, check_brenk
            _sa = calc_sa_score(mol_no_h)
            _qed = calc_qed(mol_no_h)
            _pains_ok, _pains_alerts = check_pains(mol_no_h)
            _brenk_ok, _brenk_alerts = check_brenk(mol_no_h)
            _sa_s = f"{_sa:.2f}" if not math.isnan(_sa) else "N/A"
            _qed_s = f"{_qed:.3f}" if not math.isnan(_qed) else "N/A"
            print(f"  SA={_sa_s}, QED={_qed_s}, "
                  f"PAINS={'OK' if _pains_ok else 'NG'}, "
                  f"BRENK={'OK' if _brenk_ok else 'NG'}")
        except ImportError:
            _sa = float("nan")
            _qed = float("nan")
            _pains_ok, _pains_alerts = True, []
            _brenk_ok, _brenk_alerts = True, []

        # ── ⑨ SDF 保存 ──
        sdf_path = str(OUT_DIR / f"input_{name}.sdf")
        props = {
            "AnchorA":    anchor_a["name"],
            "Scaffold":   scaffold["name"],
            "AnchorB":    anchor_b["name"],
            "Res1":       f"{resname1}{resnum1}",
            "Res2":       f"{resname2}{resnum2}",
            "TargetDist": f"{target_dist:.2f}",
            "ActualDist": f"{actual_dist:.2f}",
            "ConformanceRate": f"{conf_score['conformance_rate']:.4f}",
            "RotatableBetween": str(rot_between),
            **dl,
            "SA_Score":   f"{_sa:.2f}" if not math.isnan(_sa) else "",
            "QED":        f"{_qed:.4f}" if not math.isnan(_qed) else "",
            "PAINS_OK":   str(_pains_ok),
            "BRENK_OK":   str(_brenk_ok),
        }
        save_sdf(mol_h, best_cid, name, props, sdf_path)

        # ── ⑩ smina ドッキング ──
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
            "name":         name,
            "smiles":       canon,
            "res1":         f"{resname1}{resnum1}",
            "res2":         f"{resname2}{resnum2}",
            "anchor_a":     anchor_a["name"],
            "scaffold":     scaffold["name"],
            "anchor_b":     anchor_b["name"],
            "target_dist":  round(target_dist, 3),
            "actual_dist":  round(actual_dist, 3),
            "dock_score":   dock_score,
            "sdf_input":    sdf_path,
            "sdf_docked":   docked_sdf,
            # コンフォメーション適合性
            "conformance_rate":  conf_score["conformance_rate"],
            "median_dist":       conf_score["median_dist"],
            "dist_std":          conf_score["dist_std"],
            "rotatable_between": rot_between,
            # Drug-likeness
            **dl,
            "SA_Score":      round(_sa, 2) if not math.isnan(_sa) else None,
            "QED":           round(_qed, 4) if not math.isnan(_qed) else None,
            "PAINS_OK":      _pains_ok,
            "BRENK_OK":      _brenk_ok,
            "PAINS_Alerts":  "; ".join(_pains_alerts),
            "BRENK_Alerts":  "; ".join(_brenk_alerts),
        })

    # ── ⑪ 結果をソート・保存 ──
    # 第1ソート: conformance_rate 降順, 第2: rotatable_between 昇順, 第3: dock_score 昇順
    results.sort(key=lambda r: (
        -r.get("conformance_rate", 0.0),
        r.get("rotatable_between", 99),
        r["dock_score"] if r["dock_score"] is not None else 0.0,
    ))

    csv_path = OUT_DIR / f"rigid_{res1_str}_{res2_str}_results.csv"
    if results:
        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=results[0].keys())
            writer.writeheader()
            writer.writerows(results)
        print(f"\n結果CSV: {csv_path}")

    # ── ⑫ 可視化 ──
    if results:
        _plot_rigid_results(results, res1_str, res2_str, target_dist)

    print(f"\n{'='*60}")
    print(f"生成完了: {len(results)} 分子 (剛直スキャフォールド)")
    if results:
        best = results[0]
        print(f"最良適合率: {best['conformance_rate']:.1%}  ({best['name']})")
        if best["dock_score"] is not None:
            print(f"ドッキングスコア: {best['dock_score']:.3f} kcal/mol")
    print(f"{'='*60}\n")

    return results


# ══════════════════════════════════════════
# 可視化
# ══════════════════════════════════════════

def _plot_rigid_results(
    results: list[dict], res1: str, res2: str, target_dist: float,
):
    """剛直スキャフォールド結果の可視化"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    fig.suptitle(
        f"剛直スキャフォールドブリッジ: {res1} – {res2}  (目標 {target_dist:.1f} Å)",
        fontsize=13,
    )

    names = [r["name"].split("_")[-3] + "\n" + r["name"].split("_")[-1]
             if len(r["name"].split("_")) >= 3
             else r["name"][-20:]
             for r in results]

    # ── 適合率バー ──
    ax = axes[0]
    rates = [r.get("conformance_rate", 0) for r in results]
    colors = ["#2ca02c" if r >= 0.5 else "#ff7f0e" if r >= 0.2 else "#d62728"
              for r in rates]
    bars = ax.barh(range(len(names)), rates, color=colors)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("Conformance Rate")
    ax.set_title("ファーマコフォア適合率")
    ax.set_xlim(0, 1)
    ax.invert_yaxis()

    # ── 回転可能結合数 ──
    ax = axes[1]
    rot = [r.get("rotatable_between", 0) for r in results]
    colors_rot = ["#2ca02c" if v <= 2 else "#ff7f0e" if v <= 4 else "#d62728"
                  for v in rot]
    ax.barh(range(len(names)), rot, color=colors_rot)
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.set_xlabel("Rotatable Bonds (Between Pharmacophores)")
    ax.set_title("ファーマコフォア間回転可能結合")
    ax.invert_yaxis()

    # ── ドッキングスコア ──
    ax = axes[2]
    scores = [r["dock_score"] if r.get("dock_score") is not None else 0.0
              for r in results]
    has_dock = any(s != 0.0 for s in scores)
    if has_dock:
        colors_dock = ["#1f77b4" if s < -5 else "#ff7f0e" for s in scores]
        ax.barh(range(len(names)), scores, color=colors_dock)
        ax.set_xlabel("Docking Score (kcal/mol)")
        ax.set_title("ドッキングスコア")
    else:
        ax.text(0.5, 0.5, "ドッキングなし", ha="center", va="center",
                transform=ax.transAxes)
        ax.set_title("ドッキングスコア (N/A)")
    ax.set_yticks(range(len(names)))
    ax.set_yticklabels(names, fontsize=7)
    ax.invert_yaxis()

    plt.tight_layout()
    out_png = OUT_DIR / f"rigid_{res1}_{res2}_analysis.png"
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  分析プロット: {out_png}")


# ══════════════════════════════════════════
# CLI
# ══════════════════════════════════════════

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="剛直スキャフォールドを用いたファーマコフォアブリッジ設計",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用例:
  python rigid_scaffold_design.py --point1 LYS48 --point2 LEU50
  python rigid_scaffold_design.py --point1 LYS46 --point2 TYR49 --no-dock
  python rigid_scaffold_design.py --list-scaffolds
""",
    )
    parser.add_argument("--point1", type=str, help="残基1 (例: LYS48)")
    parser.add_argument("--point2", type=str, help="残基2 (例: LEU50)")
    parser.add_argument("--pdb", type=str, default=PDB_PATH, help="PDB ファイルパス")
    parser.add_argument("--chain", type=str, default=PROTEIN_CHAIN, help="チェーンID")
    parser.add_argument("--n-confs", type=int, default=50, help="コンフォマー生成数")
    parser.add_argument("--tolerance", type=float, default=2.0, help="距離許容誤差 (Å)")
    parser.add_argument("--exhaustiveness", type=int, default=8, help="smina exhaustiveness")
    parser.add_argument("--no-dock", action="store_true", help="ドッキングをスキップ")
    parser.add_argument("--list-scaffolds", action="store_true",
                        help="利用可能なスキャフォールドを表示")
    args = parser.parse_args()

    if args.list_scaffolds:
        _validate_scaffold_distances()
        print(f"\n{'='*60}")
        print(f"  利用可能な剛直スキャフォールド ({len(RIGID_SCAFFOLDS)} 種)")
        print(f"{'='*60}")
        print(f"  {'名前':25s} {'距離':>6s} {'回転':>4s} {'タイプ':10s} 説明")
        print(f"  {'-'*25} {'-'*6} {'-'*4} {'-'*10} {'-'*20}")
        for sc in sorted(RIGID_SCAFFOLDS, key=lambda s: s["dist_r1r2"]):
            print(f"  {sc['name']:25s} {sc['dist_r1r2']:5.1f}Å "
                  f"{sc['n_rotatable']:3d}  {sc['ring_type']:10s} {sc['desc']}")
        print(f"{'='*60}\n")
        return

    if not args.point1 or not args.point2:
        parser.error("--point1 と --point2 を指定してください")

    results = generate_rigid_bridges(
        args.point1, args.point2,
        pdb_path=args.pdb,
        chain_id=args.chain,
        dock=not args.no_dock,
        n_confs=args.n_confs,
        tolerance=args.tolerance,
        exhaustiveness=args.exhaustiveness,
    )

    # サマリー表示
    if results:
        print(f"\n{'─'*80}")
        print(f"  サマリー (上位10件)")
        print(f"{'─'*80}")
        print(f"  {'#':>3} {'名前':40s} {'適合率':>6s} {'回転':>4s} "
              f"{'スコア':>7s} {'SA':>5s} {'QED':>5s}")
        for i, r in enumerate(results[:10]):
            score_s = f"{r['dock_score']:.1f}" if r.get("dock_score") else "N/A"
            sa_s = f"{r['SA_Score']:.1f}" if r.get("SA_Score") is not None else "N/A"
            qed_s = f"{r['QED']:.3f}" if r.get("QED") is not None else "N/A"
            print(f"  {i+1:3d} {r['name']:40s} "
                  f"{r['conformance_rate']:5.1%} {r['rotatable_between']:4d} "
                  f"{score_s:>7s} {sa_s:>5s} {qed_s:>5s}")
        print(f"{'─'*80}\n")


if __name__ == "__main__":
    main()
