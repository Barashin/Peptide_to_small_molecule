#!/usr/bin/env python3
"""
compare_cyclic_peptide.py
=========================
環状ペプチド cyclo-(GEVDGWATPD) と設計済み低分子候補との結合親和性比較。

手順:
  ① 配列から head-to-tail 環状ペプチド SMILES を自動生成
  ② RDKit ETKDG + MMFF で 3D 配座を生成し SDF 保存
  ③ smina --score_only : Chain B ネイティブポーズのスコアリング
  ④ smina docking     : 環状ペプチドを de novo ドッキング
  ⑤ 既存低分子ドッキング結果 (CSV) を集約
  ⑥ 比較グラフ生成 (PRODIGY ΔG 参照線付き)

使い方:
  # ネイティブポーズ scoring + 環状ペプチドドッキング + 比較グラフ
  venv/bin/python compare_cyclic_peptide.py Protein_Peptide.pdb

  # ドッキングをスキップして既存結果と比較グラフのみ生成
  venv/bin/python compare_cyclic_peptide.py --no-dock

  # 配列を変えて比較 (デフォルト: GEVDGWATPD)
  venv/bin/python compare_cyclic_peptide.py --sequence ACDEFGHIKLM
"""

import os
import re
import sys
import csv
import json
import subprocess
import argparse
from pathlib import Path

import numpy as np

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

# 既存インフラ
import dock_with_smina as _ds


# ──────────────────────────────────────────────────────────
# 定数
# ──────────────────────────────────────────────────────────

SMINA_PATH = os.path.join(os.path.dirname(__file__), "smina.osx.12")

# ペプチド配列 (Chain B)
DEFAULT_SEQUENCE = "GEVDGWATPD"

# PRODIGY で計算済みの ΔG (kcal/mol) — 線形ペプチド複合体 (ネイティブポーズ)
PRODIGY_DG = -14.7

# 各アミノ酸の Cα 側鎖 SMILES (バックボーン N, Cα, CO は別途構築)
_AA_SIDE: dict[str, str | None] = {
    "G": "",              # Glycine: 側鎖なし
    "A": "C",             # Alanine
    "V": "C(C)C",         # Valine
    "L": "CC(C)C",        # Leucine
    "I": "C(C)CC",        # Isoleucine
    "F": "Cc1ccccc1",     # Phenylalanine
    "W": "Cc1c[nH]c2ccccc12",  # Tryptophan
    "M": "CCSC",          # Methionine
    "S": "CO",            # Serine
    "T": "C(O)C",         # Threonine
    "C": "CS",            # Cysteine
    "Y": "Cc1ccc(O)cc1",  # Tyrosine
    "H": "Cc1cnc[nH]1",   # Histidine
    "D": "CC(=O)O",       # Aspartate
    "E": "CCC(=O)O",      # Glutamate
    "N": "CC(=O)N",       # Asparagine
    "Q": "CCC(=O)N",      # Glutamine
    "K": "CCCCN",         # Lysine
    "R": "CCCNC(=N)N",    # Arginine
    "P": None,            # Proline: 特殊 (ピロリジン環)
}


# ──────────────────────────────────────────────────────────
# 1. 環状ペプチド SMILES 生成
# ──────────────────────────────────────────────────────────

def build_cyclic_peptide_smiles(sequence: str) -> str:
    """
    head-to-tail 環状ペプチド SMILES を生成する。

    バックボーンパターン (各残基):
      通常残基: N-Cα(side)-C(=O)
      Gly    : N-C-C(=O)  (側鎖なし)
      Pro    : N{ring}CCC{ring}C(=O)  (ピロリジン環)

    環化:
      最初の残基の N に ring closure %11 を開き、
      最後の残基の CO で C(=O)%11 として閉じる。

    立体化学は省略 (比較には不要)。
    """
    seq = sequence.upper()
    n   = len(seq)
    if n < 2:
        raise ValueError("環状ペプチドには2残基以上必要です")

    # head-to-tail リング番号 (%11)
    HEAD_RING = "%11"
    # Pro ピロリジン環番号 (3, 4, 5 ... 順に割り当て; 1,2 は Trp/His の芳香環)
    pro_ring_counter = [3]

    def _residue_smiles(aa: str, pos: str) -> str:
        """'first' / 'middle' / 'last' に応じたバックボーンフラグメントを返す"""
        side = _AA_SIDE.get(aa, "")

        if aa == "P":
            r = pro_ring_counter[0]
            pro_ring_counter[0] += 1
            # Pro: バックボーン N が ピロリジン環の一部
            # N{r}CCCC{r}C(=O) — N{r}開環, CCCC{r}=Cδ,Cγ,Cβ,Cα(閉環), CO
            if pos == "first":
                # head-to-tail 開環が N にある場合: Pro は先頭にできない制約
                raise ValueError(
                    "Pro が先頭の環状ペプチドは SMILES 生成が複雑なため "
                    "--sequence で別の残基を先頭にしてください"
                )
            core = f"N{r}CCCC{r}C(=O)"
            return core

        # 通常残基
        if side:
            ca_part = f"C({side})"
        else:
            ca_part = "C"  # Gly

        if pos == "first":
            return f"N{HEAD_RING}{ca_part}C(=O)"
        elif pos == "last":
            return f"N{ca_part}C(=O){HEAD_RING}"
        else:
            return f"N{ca_part}C(=O)"

    parts = []
    for i, aa in enumerate(seq):
        if i == 0:
            pos = "first"
        elif i == n - 1:
            pos = "last"
        else:
            pos = "middle"
        parts.append(_residue_smiles(aa, pos))

    return "".join(parts)


def smiles_to_mol(smiles: str) -> Chem.Mol | None:
    """SMILES → RDKit Mol (バリデーション付き)"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  [警告] SMILES のパース失敗: {smiles[:80]}...")
    return mol


# ──────────────────────────────────────────────────────────
# Ligand Efficiency (LE) 計算
# ──────────────────────────────────────────────────────────

def calc_hac(smiles: str) -> int | None:
    """SMILES から重原子数 (Heavy Atom Count) を返す"""
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms() if mol is not None else None


def calc_le(score: float | None, hac: int | None) -> float | None:
    """
    Ligand Efficiency = |ΔG (kcal/mol)| / HAC

    LE > 0.3  : 効率的な結合剤 (低分子医薬の目安)
    LE > 0.2  : 許容範囲
    LE < 0.2  : 非効率 (大きすぎ or 弱すぎ)
    """
    if score is None or hac is None or hac == 0:
        return None
    return round(abs(score) / hac, 4)


# ──────────────────────────────────────────────────────────
# 2. 3D 配座生成
# ──────────────────────────────────────────────────────────

def generate_conformer(mol: Chem.Mol,
                       n_confs: int = 300,
                       seed: int = 42) -> Chem.Mol | None:
    """
    ETKDG で大環状の 3D 配座を生成し、MMFF 最適化後に最低エネルギー配座を返す。
    ETKDG が全失敗した場合は DG (classic) にフォールバック。
    """
    mol_h = Chem.AddHs(mol)

    # ETKDG パラメータ (大環状対応)
    params = AllChem.ETKDGv3()
    params.randomSeed         = seed
    params.useSmallRingTorsions = True
    params.useMacrocycleTorsions = True  # RDKit >= 2021 で大環状をサポート
    params.numThreads         = 0  # マルチスレッド

    confs = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params)
    n_generated = len(confs)
    print(f"  ETKDG 生成配座数: {n_generated} / {n_confs}")

    if n_generated == 0:
        # フォールバック: classic DG
        print("  ETKDG 失敗 → classic DG にフォールバック")
        params_dg = AllChem.EmbedParameters()
        params_dg.randomSeed = seed
        confs = AllChem.EmbedMultipleConfs(mol_h, numConfs=n_confs, params=params_dg)
        n_generated = len(confs)
        if n_generated == 0:
            print("  [エラー] 3D 配座生成に完全失敗しました")
            return None

    # MMFF 力場最適化
    results = AllChem.MMFFOptimizeMoleculeConfs(mol_h, mmffVariant="MMFF94s",
                                                maxIters=2000, numThreads=0)

    # 最低 MMFF エネルギーの配座を選択
    best_conf, best_energy = -1, float("inf")
    for conf_id, (converged, energy) in enumerate(results):
        if energy < best_energy:
            best_energy = energy
            best_conf = conf_id

    if best_conf == -1:
        best_conf = 0

    # 選択された配座だけを保持
    mol_out = Chem.RWMol(mol_h)
    confs_to_remove = [i for i in range(mol_h.GetNumConformers()) if i != best_conf]
    for c in sorted(confs_to_remove, reverse=True):
        mol_out.RemoveConformer(c)

    mol_final = Chem.RemoveHs(mol_out)
    print(f"  最良配座: #{best_conf}, MMFF エネルギー = {best_energy:.2f}")
    return mol_final


def write_mol_to_sdf(mol: Chem.Mol, sdf_path: str, name: str = "cyclic_peptide"):
    """Mol を SDF ファイルに書き出す"""
    mol.SetProp("_Name", name)
    smiles = Chem.MolToSmiles(mol)
    mol.SetProp("SMILES", smiles)
    w = Chem.SDWriter(sdf_path)
    w.write(mol)
    w.close()


# ──────────────────────────────────────────────────────────
# 3. ネイティブポーズ scoring (smina --score_only)
# ──────────────────────────────────────────────────────────

def score_native_pose(peptide_pdb: str,
                      receptor_pdb: str,
                      log_path: str) -> float | None:
    """
    Chain B のネイティブ結合ポーズを smina --score_only でスコアリング。
    peptide_pdb には結合状態の Chain B PDB を渡す。
    """
    print(f"\n  [score_only] ネイティブポーズをスコアリング中...")
    print(f"    受容体 : {receptor_pdb}")
    print(f"    リガンド: {peptide_pdb}")

    cmd = [
        SMINA_PATH,
        "--receptor",   receptor_pdb,
        "--ligand",     peptide_pdb,
        "--score_only",
        "--log",        log_path,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    stdout = result.stdout + result.stderr

    affinity_m = re.search(r"Affinity:\s*([-\d.]+)", stdout)
    if affinity_m:
        score = float(affinity_m.group(1))
        print(f"    → ネイティブポーズ score_only: {score:.3f} kcal/mol")
        return score

    print("  [警告] ネイティブポーズスコアリング失敗")
    print(stdout[-1000:])
    return None


# ──────────────────────────────────────────────────────────
# 4. smina ドッキング
# ──────────────────────────────────────────────────────────

def dock_single_mol(sdf_path: str,
                    receptor_pdb: str,
                    ref_pdb: str,
                    output_sdf: str,
                    log_path: str,
                    exhaustiveness: int = 16,
                    num_modes: int = 9) -> float | None:
    """
    単一分子 SDF をドッキングし、ベストスコアを返す。
    環状ペプチドは自由度が高いため exhaustiveness=16 以上を推奨。
    """
    print(f"\n  [docking] {Path(sdf_path).name} をドッキング中...")
    stdout = _ds.run_smina(
        receptor_pdb, sdf_path, ref_pdb,
        output_sdf, log_path,
        exhaustiveness=exhaustiveness,
        num_modes=num_modes,
        autobox_add=4.0,
    )
    scores = _ds.parse_smina_scores(log_path, stdout)
    if scores:
        best = min(scores)
        print(f"    → Best docking score: {best:.3f} kcal/mol")
        return best
    print("  [警告] ドッキングスコアが取得できませんでした")
    return None


# ──────────────────────────────────────────────────────────
# 5. 既存低分子ドッキング結果の集約
# ──────────────────────────────────────────────────────────

def load_peptidomimetics_results(results_dir: str,
                                 max_n: int = 15) -> list[dict]:
    """results/docking/docking_results.csv から読み込む"""
    csv_path = os.path.join(results_dir, "docking", "docking_results.csv")
    if not os.path.exists(csv_path):
        print(f"  [スキップ] {csv_path} が見つかりません")
        return []

    records = []
    with open(csv_path) as f:
        for row in csv.DictReader(f):
            score_str = row.get("best_score_kcal/mol", "N/A")
            if score_str == "N/A":
                continue
            try:
                score  = float(score_str)
                smiles = row.get("smiles", "")
                hac    = calc_hac(smiles)
                records.append({
                    "name"  : row["name"],
                    "score" : score,
                    "hac"   : hac,
                    "le"    : calc_le(score, hac),
                    "smiles": smiles,
                    "source": "ペプチドミメティクス (STEP3)",
                    "color" : "#3498db",
                })
            except ValueError:
                pass

    records.sort(key=lambda r: r["score"])
    print(f"  ペプチドミメティクス: {len(records)} 件 読み込み")
    return records[:max_n]


def load_bridge_results(results_dir: str,
                        max_n: int = 15) -> list[dict]:
    """results/bridge/bridge_*_results.csv から読み込む"""
    bridge_dir = os.path.join(results_dir, "bridge")
    if not os.path.isdir(bridge_dir):
        print(f"  [スキップ] {bridge_dir} が見つかりません")
        return []

    all_records = []
    for csv_path in Path(bridge_dir).glob("bridge_*_results.csv"):
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                score_str = row.get("dock_score", "").strip()
                if not score_str:
                    continue
                try:
                    score  = float(score_str)
                    smiles = row.get("smiles", "")
                    hac    = calc_hac(smiles)
                    all_records.append({
                        "name"  : row["name"],
                        "score" : score,
                        "hac"   : hac,
                        "le"    : calc_le(score, hac),
                        "smiles": smiles,
                        "source": "ファーマコフォアブリッジ (STEP7)",
                        "color" : "#e67e22",
                    })
                except ValueError:
                    pass

    all_records.sort(key=lambda r: r["score"])
    print(f"  ファーマコフォアブリッジ: {len(all_records)} 件 読み込み")
    return all_records[:max_n]


# ──────────────────────────────────────────────────────────
# 6. 比較グラフ生成
# ──────────────────────────────────────────────────────────

def plot_comparison(peptide_entries: list[dict],
                    ligand_entries: list[dict],
                    prodigy_dg: float,
                    output_path: str):
    """
    環状ペプチド / ネイティブポーズ / 低分子候補の結合親和性比較グラフ。

    Parameters
    ----------
    peptide_entries : [{"label", "score", "color"}, ...]  # ペプチド系
    ligand_entries  : [{"name", "score", "source", "color"}, ...]
    prodigy_dg      : PRODIGY で計算した ΔG (kcal/mol)
    output_path     : 出力 PNG パス
    """
    import plot_utils  # noqa: F401  (日本語フォント設定)
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    # ─── データ準備 ───
    all_items = []

    # ペプチド系 (上部に配置)
    for p in peptide_entries:
        all_items.append({
            "label": p["label"],
            "score": p["score"],
            "color": p["color"],
            "hatch": p.get("hatch", ""),
        })

    # 区切り線用ダミー (プロット上は除外して処理)
    SEPARATOR = "─────────"

    # 低分子候補 (スコア昇順)
    sorted_ligands = sorted(ligand_entries, key=lambda r: r["score"])
    for r in sorted_ligands:
        all_items.append({
            "label": r["name"][:35],
            "score": r["score"],
            "color": r["color"],
            "hatch": "",
        })

    if not all_items:
        print("  [警告] 比較するデータがありません")
        return

    labels = [it["label"] for it in all_items]
    scores = [it["score"] for it in all_items]
    colors = [it["color"] for it in all_items]

    n = len(all_items)
    fig_h = max(6, n * 0.55 + 3)
    fig, ax = plt.subplots(figsize=(12, fig_h))

    y_pos = np.arange(n)
    abs_scores = [abs(s) for s in scores]

    bars = ax.barh(y_pos, abs_scores, color=colors, edgecolor="white",
                   linewidth=0.8, height=0.75)

    # スコア値ラベル
    for bar, s in zip(bars, scores):
        ax.text(bar.get_width() + 0.08, bar.get_y() + bar.get_height() / 2,
                f"{s:.2f}", va="center", ha="left", fontsize=9)

    # ─── PRODIGY 参照線 ───
    ax.axvline(abs(prodigy_dg), color="#c0392b", linewidth=1.5,
               linestyle="--", alpha=0.8,
               label=f"PRODIGY ΔG = {prodigy_dg:.1f} kcal/mol\n(線形ペプチド, 別スケール参考)")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("|Affinity| (kcal/mol)  — 右ほど結合が強い", fontsize=11)
    ax.set_title("結合親和性比較: 環状ペプチド vs 設計低分子", fontsize=14, fontweight="bold")

    # ─── 凡例 ───
    legend_handles = [
        mpatches.Patch(color="#2ecc71", label="環状ペプチド (smina dock)"),
        mpatches.Patch(color="#27ae60", label="線形ペプチド (smina score_only)"),
        mpatches.Patch(color="#3498db", label="ペプチドミメティクス (STEP3)"),
        mpatches.Patch(color="#e67e22", label="ファーマコフォアブリッジ (STEP7)"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=9)

    # ─── ペプチド vs 低分子の区切り線 ───
    if peptide_entries:
        n_pep = len(peptide_entries)
        ax.axhline(n_pep - 0.5, color="#7f8c8d", linewidth=1.0, linestyle=":")

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"\n  比較グラフ 保存: {output_path}")


# ──────────────────────────────────────────────────────────
# 7. Ligand Efficiency 比較グラフ
# ──────────────────────────────────────────────────────────

def plot_le_comparison(peptide_entries: list[dict],
                       ligand_entries: list[dict],
                       output_path: str):
    """
    Ligand Efficiency (LE = |ΔG| / HAC) 比較グラフ。

    LE で比較することで分子サイズのバイアスを除去し、
    ペプチド (大分子) と低分子を公平に評価できる。

    参考ライン:
      LE = 0.3  (低分子医薬品として効率的な結合の目安)
      LE = 0.2  (許容下限)
    """
    import plot_utils  # noqa: F401
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches

    all_items = []

    for p in peptide_entries:
        if p.get("le") is not None:
            all_items.append({
                "label": p["label"],
                "le"   : p["le"],
                "hac"  : p.get("hac", "?"),
                "score": p["score"],
                "color": p["color"],
            })

    sorted_ligands = sorted(
        [r for r in ligand_entries if r.get("le") is not None],
        key=lambda r: r["le"],
        reverse=True,   # LE 高い順 (上が優秀)
    )
    for r in sorted_ligands:
        all_items.append({
            "label": r["name"][:35],
            "le"   : r["le"],
            "hac"  : r.get("hac", "?"),
            "score": r["score"],
            "color": r["color"],
        })

    if not all_items:
        print("  [警告] LE を計算できるデータがありません")
        return

    labels = [it["label"] for it in all_items]
    les    = [it["le"]    for it in all_items]
    colors = [it["color"] for it in all_items]
    hacs   = [it["hac"]   for it in all_items]

    n     = len(all_items)
    fig_h = max(6, n * 0.55 + 3)
    fig, ax = plt.subplots(figsize=(12, fig_h))

    y_pos = np.arange(n)
    bars  = ax.barh(y_pos, les, color=colors, edgecolor="white",
                    linewidth=0.8, height=0.75)

    # LE 値 + HAC ラベル
    for bar, le, hac in zip(bars, les, hacs):
        ax.text(bar.get_width() + 0.002,
                bar.get_y() + bar.get_height() / 2,
                f"LE={le:.3f}  (HAC={hac})",
                va="center", ha="left", fontsize=8.5)

    # ─── 参照ライン ───
    ax.axvline(0.3, color="#27ae60", linewidth=1.5, linestyle="--",
               label="LE = 0.3 (効率的な結合の目安)")
    ax.axvline(0.2, color="#e67e22", linewidth=1.2, linestyle=":",
               label="LE = 0.2 (許容下限)")

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Ligand Efficiency  LE = |ΔG| / HAC  (kcal/mol/atom)",
                  fontsize=11)
    ax.set_title("Ligand Efficiency 比較\n"
                 "─ サイズバイアスを除いた公平な結合効率評価 ─",
                 fontsize=13, fontweight="bold")

    # ─── 凡例 ───
    legend_handles = [
        mpatches.Patch(color="#2ecc71", label="環状ペプチド (smina dock)"),
        mpatches.Patch(color="#27ae60", label="線形ペプチド (smina score_only)"),
        mpatches.Patch(color="#3498db", label="ペプチドミメティクス (STEP3)"),
        mpatches.Patch(color="#e67e22", label="ファーマコフォアブリッジ (STEP7)"),
    ]
    ax.legend(handles=legend_handles, loc="lower right", fontsize=9)

    # ─── ペプチド vs 低分子の区切り線 ───
    n_pep = sum(1 for p in peptide_entries if p.get("le") is not None)
    if n_pep > 0:
        ax.axhline(n_pep - 0.5, color="#7f8c8d", linewidth=1.0, linestyle=":")

    # ─── 注釈 ───
    ax.text(0.01, -0.06,
            "※ LE = Ligand Efficiency = |smina Affinity| / 重原子数。"
            "分子サイズの影響を除いて結合効率を評価する指標。",
            transform=ax.transAxes, fontsize=8, color="#555555")

    plt.tight_layout(rect=[0, 0.04, 1, 1])
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  LE 比較グラフ 保存: {output_path}")


# ──────────────────────────────────────────────────────────
# 8. 結果サマリー表示
# ──────────────────────────────────────────────────────────

def print_summary(peptide_entries: list[dict],
                  ligand_entries: list[dict],
                  prodigy_dg: float):
    print()
    print("=" * 90)
    print("  結合親和性比較サマリー")
    print(f"  {'名称':<42} {'Score':>9}  {'HAC':>5}  {'LE':>7}  種別")
    print("=" * 90)

    print(f"\n  【参照: ペプチド系 (smina)】")
    for p in peptide_entries:
        score_str = f"{p['score']:.3f}" if p.get("score") is not None else "N/A"
        hac_str   = str(p.get("hac", "?"))
        le_str    = f"{p['le']:.4f}" if p.get("le") is not None else "N/A"
        print(f"  ★ {p['label']:<40} {score_str:>9}  {hac_str:>5}  {le_str:>7}")

    print(f"\n  【参照: PRODIGY ΔG (線形ペプチド複合体, 別スケール)】")
    print(f"  {'PRODIGY':<42} {prodigy_dg:>9.1f}  {'74':>5}  {'--':>7}  protein-peptide 統計モデル")
    print(f"    ※ smina Vina スコアとは定義・スケールが異なるため LE は非表示")

    if ligand_entries:
        print(f"\n  【設計低分子 (smina docking, スコア上位15件)】")
        ranked = sorted(ligand_entries, key=lambda r: r["score"])
        for rank, r in enumerate(ranked[:15], 1):
            hac_str = str(r.get("hac", "?"))
            le_str  = f"{r['le']:.4f}" if r.get("le") is not None else "N/A"
            print(f"  {rank:>2}. {r['name']:<40} {r['score']:>9.3f}  {hac_str:>5}  {le_str:>7}  {r['source']}")

    print()
    print("  LE = Ligand Efficiency = |Score| / HAC  (kcal/mol/重原子)")
    print("  LE > 0.3: 効率的な結合剤 (低分子医薬の目安)")
    print("  LE > 0.2: 許容範囲  |  LE < 0.2: 非効率 (大きすぎ or 弱すぎ)")
    print("=" * 90)


# ──────────────────────────────────────────────────────────
# 8. 分子プロパティ表示
# ──────────────────────────────────────────────────────────

def print_mol_properties(mol: Chem.Mol, name: str):
    """環状ペプチドの物理化学的性質を表示"""
    mw   = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd  = rdMolDescriptors.CalcNumHBD(mol)
    hba  = rdMolDescriptors.CalcNumHBA(mol)
    rot  = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)
    print(f"\n  [{name}] 物理化学的性質:")
    print(f"    MW = {mw:.1f} Da  |  LogP = {logp:.2f}  |  "
          f"HBD = {hbd}  |  HBA = {hba}  |  RotBonds = {rot}  |  Rings = {rings}")


# ──────────────────────────────────────────────────────────
# メイン
# ──────────────────────────────────────────────────────────

def main():
    p = argparse.ArgumentParser(
        description="環状ペプチドと設計低分子の結合親和性比較")
    p.add_argument("pdb", nargs="?", default="Protein_Peptide.pdb",
                   help="入力 PDB ファイル (デフォルト: Protein_Peptide.pdb)")
    p.add_argument("--sequence",   default=DEFAULT_SEQUENCE,
                   help=f"ペプチド配列 1文字コード (デフォルト: {DEFAULT_SEQUENCE})")
    p.add_argument("--no-dock",    action="store_true",
                   help="smina ドッキングをスキップ (グラフ生成のみ)")
    p.add_argument("--n-confs",    type=int, default=300,
                   help="3D 配座生成数 (デフォルト: 300)")
    p.add_argument("--exhaustiveness", type=int, default=16,
                   help="smina 探索の徹底度 (デフォルト: 16)")
    p.add_argument("--results-dir", default="results",
                   help="既存結果ディレクトリ (デフォルト: results)")
    p.add_argument("--output-dir", default="results/comparison",
                   help="出力先ディレクトリ (デフォルト: results/comparison)")
    p.add_argument("--protein-chain", default="A",
                   help="タンパク質チェーン ID (デフォルト: A)")
    p.add_argument("--peptide-chain", default="B",
                   help="ペプチドチェーン ID (デフォルト: B)")
    args = p.parse_args()

    out_dir = args.output_dir
    os.makedirs(out_dir, exist_ok=True)
    res_dir = args.results_dir

    print("\n" + "=" * 60)
    print("  環状ペプチド vs 低分子 結合親和性比較ツール")
    print("=" * 60)
    print(f"  配列    : {args.sequence}")
    print(f"  PDB     : {args.pdb}")
    print(f"  出力先  : {out_dir}")

    # ─── 受容体・参照ペプチドの準備 ───
    receptor_pdb = os.path.join(res_dir, "docking", "receptor.pdb")
    peptide_ref  = os.path.join(res_dir, "docking", "peptide_ref.pdb")

    if not os.path.exists(receptor_pdb) or not os.path.exists(peptide_ref):
        print("\n  受容体/参照ペプチドが未生成 → PDB からチェーン抽出")
        dock_dir = os.path.join(res_dir, "docking")
        os.makedirs(dock_dir, exist_ok=True)
        _ds.extract_chain(args.pdb, args.protein_chain, receptor_pdb)
        _ds.extract_chain(args.pdb, args.peptide_chain, peptide_ref)

    # ─── 環状ペプチド SMILES 生成 ───
    print(f"\n  ① 環状ペプチド SMILES 生成 (cyclo-{args.sequence})")
    try:
        cyclic_smiles = build_cyclic_peptide_smiles(args.sequence)
    except ValueError as e:
        print(f"  [エラー] {e}")
        sys.exit(1)

    print(f"  SMILES: {cyclic_smiles[:100]}{'...' if len(cyclic_smiles) > 100 else ''}")

    # SMILES バリデーション
    cyclic_mol = smiles_to_mol(cyclic_smiles)
    if cyclic_mol is None:
        print("  [エラー] SMILES のパースに失敗しました。配列を確認してください。")
        sys.exit(1)

    print_mol_properties(cyclic_mol, f"cyclo-{args.sequence}")

    # ─── 3D 配座生成 ───
    cyclic_sdf  = os.path.join(out_dir, f"cyclic_peptide_{args.sequence}.sdf")
    cyclic_name = f"cyclo_{args.sequence}"

    if not args.no_dock:
        print(f"\n  ② 3D 配座生成 (ETKDG + MMFF, n_confs={args.n_confs})")
        mol_3d = generate_conformer(cyclic_mol, n_confs=args.n_confs)
        if mol_3d is None:
            print("  [エラー] 3D 配座生成に失敗しました")
            sys.exit(1)
        write_mol_to_sdf(mol_3d, cyclic_sdf, name=cyclic_name)
        print(f"  SDF 保存: {cyclic_sdf}")

    # ─── ネイティブポーズ scoring ───
    native_score = None
    if not args.no_dock:
        print(f"\n  ③ 線形ペプチド ネイティブポーズ score_only")
        native_log = os.path.join(out_dir, "native_peptide_score.log")
        native_score = score_native_pose(peptide_ref, receptor_pdb, native_log)

    # ─── 環状ペプチド ドッキング ───
    cyclic_score = None
    if not args.no_dock:
        print(f"\n  ④ 環状ペプチド ドッキング (exhaustiveness={args.exhaustiveness})")
        cyclic_docked_sdf = os.path.join(out_dir, f"docked_{cyclic_name}.sdf")
        cyclic_log        = os.path.join(out_dir, f"docked_{cyclic_name}.log")
        cyclic_score = dock_single_mol(
            cyclic_sdf, receptor_pdb, peptide_ref,
            cyclic_docked_sdf, cyclic_log,
            exhaustiveness=args.exhaustiveness,
        )
    else:
        # --no-dock の場合: 以前の結果が残っていれば読み込む
        prev_log = os.path.join(out_dir, f"docked_{cyclic_name}.log")
        prev_sdf = os.path.join(out_dir, f"docked_{cyclic_name}.sdf")
        if os.path.exists(prev_log):
            with open(prev_log) as f:
                stdout = f.read()
            scores = _ds.parse_smina_scores(prev_log, stdout)
            if scores:
                cyclic_score = min(scores)
                print(f"  [前回結果] 環状ペプチド docking score: {cyclic_score:.3f}")
        native_log = os.path.join(out_dir, "native_peptide_score.log")
        if os.path.exists(native_log):
            with open(native_log) as f:
                stdout = f.read()
            m = re.search(r"Affinity:\s*([-\d.]+)", stdout)
            if m:
                native_score = float(m.group(1))
                print(f"  [前回結果] 線形ペプチド score_only: {native_score:.3f}")

    # ─── 既存低分子結果の集約 ───
    print(f"\n  ⑤ 既存低分子ドッキング結果の集約")
    peptidomimetics = load_peptidomimetics_results(res_dir, max_n=20)
    bridge_mols     = load_bridge_results(res_dir, max_n=20)
    all_ligands     = peptidomimetics + bridge_mols

    # ─── ペプチド系エントリ (HAC・LE を付与) ───
    # 環状ペプチド HAC: RDKit で計算済み
    cyclic_hac  = cyclic_mol.GetNumAtoms()
    # 線形ペプチド HAC: 環状 +1 (head-to-tail 閉環で O が 1 個余分に失われるため)
    linear_hac  = cyclic_hac + 1

    peptide_entries = []
    if cyclic_score is not None:
        peptide_entries.append({
            "label": f"cyclo-{args.sequence} (smina dock)",
            "score": cyclic_score,
            "hac"  : cyclic_hac,
            "le"   : calc_le(cyclic_score, cyclic_hac),
            "color": "#2ecc71",
        })
    if native_score is not None:
        peptide_entries.append({
            "label": f"linear {args.sequence} (native, score_only)",
            "score": native_score,
            "hac"  : linear_hac,
            "le"   : calc_le(native_score, linear_hac),
            "color": "#27ae60",
        })

    # ─── 結果サマリー ───
    print_summary(peptide_entries, all_ligands, PRODIGY_DG)

    # ─── 比較グラフ生成 ───
    print(f"\n  ⑥ 比較グラフ生成")
    # 上位低分子に絞る (スコア上位 15 件)
    top_ligands = sorted(all_ligands, key=lambda r: r["score"])[:15]

    graph_path = os.path.join(out_dir, "binding_comparison.png")
    plot_comparison(peptide_entries, top_ligands, PRODIGY_DG, graph_path)

    le_graph_path = os.path.join(out_dir, "le_comparison.png")
    plot_le_comparison(peptide_entries, top_ligands, le_graph_path)

    # ─── JSON サマリー保存 ───
    summary = {
        "sequence"        : args.sequence,
        "prodigy_dg"      : PRODIGY_DG,
        "cyclic_peptide"  : {
            "smiles"    : cyclic_smiles,
            "hac"       : cyclic_hac,
            "dock_score": cyclic_score,
            "le"        : calc_le(cyclic_score, cyclic_hac),
        },
        "linear_peptide"  : {
            "hac"       : linear_hac,
            "score_only": native_score,
            "le"        : calc_le(native_score, linear_hac),
        },
        "top_ligands"     : [
            {
                "name"  : r["name"],
                "score" : r["score"],
                "hac"   : r.get("hac"),
                "le"    : r.get("le"),
                "source": r["source"],
            }
            for r in top_ligands
        ],
    }
    json_path = os.path.join(out_dir, "comparison_summary.json")
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    print(f"  JSON サマリー保存: {json_path}")

    print("\n  完了!")


if __name__ == "__main__":
    main()
