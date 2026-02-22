#!/usr/bin/env python3
"""
synthesizability.py
===================
合成容易性 (SA) スコア・構造アラート・薬剤らしさ (QED) の計算モジュール。

提供する機能:
  - SA Score  (合成容易性、1〜10、低い = 合成しやすい)
  - PAINS フィルタ  (pan-assay interference compounds, 480 パターン)
  - BRENK フィルタ  (構造アラート, 105 パターン)
  - QED  (Quantitative Estimate of Drug-likeness, 0〜1、高い = 良)
  - バッチエンリッチ関数  (候補 dict にスコアを追加)
  - フィルタ関数  (閾値ベースで候補を選別)
"""

import math
import sys
import os

from rdkit import Chem
from rdkit.Chem import QED as _QED_module
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


# ──────────────────────────────────────────────────────────
# SA Score (lazy-loaded singleton)
# ──────────────────────────────────────────────────────────

_SASCORER = None


def _load_sascorer():
    """RDKit Contrib/SA_Score/sascorer.py を動的にインポート"""
    try:
        from rdkit.Chem import RDConfig
        sa_dir = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if sa_dir not in sys.path:
            sys.path.insert(0, sa_dir)
        import sascorer
        return sascorer
    except Exception as e:
        print(f"  [警告] SA Score モジュールの読み込みに失敗: {e}")
        return None


def calc_sa_score(mol) -> float:
    """SA Score を計算 (1.0 = 合成容易, 10.0 = 合成困難)。失敗時は nan。"""
    global _SASCORER
    if _SASCORER is None:
        _SASCORER = _load_sascorer()
    if _SASCORER is None or mol is None:
        return float("nan")
    try:
        return round(_SASCORER.calculateScore(mol), 2)
    except Exception:
        return float("nan")


# ──────────────────────────────────────────────────────────
# PAINS / BRENK フィルタ (lazy-loaded singletons)
# ──────────────────────────────────────────────────────────

_PAINS_CAT = None
_BRENK_CAT = None


def _get_pains_catalog() -> FilterCatalog:
    global _PAINS_CAT
    if _PAINS_CAT is None:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
        _PAINS_CAT = FilterCatalog(params)
    return _PAINS_CAT


def _get_brenk_catalog() -> FilterCatalog:
    global _BRENK_CAT
    if _BRENK_CAT is None:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        _BRENK_CAT = FilterCatalog(params)
    return _BRENK_CAT


def check_pains(mol) -> tuple:
    """
    PAINS フィルタ適用。
    Returns: (is_clean: bool, alert_names: list[str])
        is_clean = True  → PAINS アラートなし (良好)
    """
    if mol is None:
        return False, ["invalid_mol"]
    cat = _get_pains_catalog()
    entry = cat.GetFirstMatch(mol)
    if entry is None:
        return True, []
    alerts = [m.GetDescription() for m in cat.GetMatches(mol)]
    return False, alerts


def check_brenk(mol) -> tuple:
    """
    BRENK 構造アラートフィルタ適用。
    Returns: (is_clean: bool, alert_names: list[str])
    """
    if mol is None:
        return False, ["invalid_mol"]
    cat = _get_brenk_catalog()
    entry = cat.GetFirstMatch(mol)
    if entry is None:
        return True, []
    alerts = [m.GetDescription() for m in cat.GetMatches(mol)]
    return False, alerts


# ──────────────────────────────────────────────────────────
# QED
# ──────────────────────────────────────────────────────────

def calc_qed(mol) -> float:
    """QED (Quantitative Estimate of Drug-likeness)。0〜1 (高い = 良好)。"""
    if mol is None:
        return float("nan")
    try:
        return round(_QED_module.qed(mol), 4)
    except Exception:
        return float("nan")


# ──────────────────────────────────────────────────────────
# バッチエンリッチ
# ──────────────────────────────────────────────────────────

def enrich_candidates(candidates: list, verbose: bool = True) -> list:
    """
    候補分子 dict のリストに合成容易性スコアを in-place で追加する。

    追加されるキー:
      SA_Score     : float  (1〜10, 低い = 合成しやすい)
      QED          : float  (0〜1, 高い = 薬剤らしい)
      PAINS_OK     : bool   (True = PAINS アラートなし)
      PAINS_Alerts : str    (セミコロン区切り)
      BRENK_OK     : bool   (True = BRENK アラートなし)
      BRENK_Alerts : str
      Synth_Pass   : bool   (SA ≤ 6.0 and PAINS_OK and BRENK_OK)

    Args:
        candidates: 各 dict に "mol" キー (rdkit.Chem.Mol) が必要。
                    "mol" がなく "smiles" がある場合はそこから生成する。
        verbose: True なら合計統計を表示
    """
    n_pains = 0
    n_brenk = 0
    n_sa_hard = 0

    for c in candidates:
        mol = c.get("mol")
        if mol is None and c.get("smiles"):
            mol = Chem.MolFromSmiles(c["smiles"])

        # SA Score
        sa = calc_sa_score(mol)
        c["SA_Score"] = sa

        # QED
        c["QED"] = calc_qed(mol)

        # PAINS
        pains_ok, pains_alerts = check_pains(mol)
        c["PAINS_OK"] = pains_ok
        c["PAINS_Alerts"] = "; ".join(pains_alerts)
        if not pains_ok:
            n_pains += 1

        # BRENK
        brenk_ok, brenk_alerts = check_brenk(mol)
        c["BRENK_OK"] = brenk_ok
        c["BRENK_Alerts"] = "; ".join(brenk_alerts)
        if not brenk_ok:
            n_brenk += 1

        # 合成容易性合格判定
        sa_ok = (not math.isnan(sa)) and sa <= 6.0
        c["Synth_Pass"] = sa_ok and pains_ok and brenk_ok
        if not sa_ok:
            n_sa_hard += 1

    if verbose and candidates:
        total = len(candidates)
        n_pass = sum(1 for c in candidates if c.get("Synth_Pass", False))
        print(f"\n  [合成容易性] {total} 分子を評価:")
        print(f"    SA Score > 6.0 (合成困難): {n_sa_hard}/{total}")
        print(f"    PAINS アラート:            {n_pains}/{total}")
        print(f"    BRENK アラート:            {n_brenk}/{total}")
        print(f"    合格 (SA≤6 & PAINS_OK & BRENK_OK): {n_pass}/{total}")

    return candidates


# ──────────────────────────────────────────────────────────
# フィルタ
# ──────────────────────────────────────────────────────────

DEFAULT_THRESHOLDS = {
    "SA_Score_max": 6.0,
    "QED_min":      0.2,
    "PAINS_OK":     True,
    "BRENK_OK":     True,
}


def filter_candidates(
    candidates: list,
    thresholds: dict = None,
    verbose: bool = True,
) -> list:
    """
    閾値ベースで候補分子をフィルタする。

    Returns:
        フィルタ通過した候補のリスト (元リストは変更しない)
    """
    th = {**DEFAULT_THRESHOLDS, **(thresholds or {})}
    passed = []
    for c in candidates:
        sa = c.get("SA_Score", float("nan"))
        qed = c.get("QED", float("nan"))

        if not math.isnan(sa) and sa > th["SA_Score_max"]:
            continue
        if not math.isnan(qed) and qed < th["QED_min"]:
            continue
        if th["PAINS_OK"] and not c.get("PAINS_OK", True):
            continue
        if th["BRENK_OK"] and not c.get("BRENK_OK", True):
            continue
        passed.append(c)

    if verbose:
        print(f"  [フィルタ] {len(candidates)} → {len(passed)} 分子 "
              f"(SA≤{th['SA_Score_max']}, QED≥{th['QED_min']}, "
              f"PAINS={'除外' if th['PAINS_OK'] else '許容'}, "
              f"BRENK={'除外' if th['BRENK_OK'] else '許容'})")

    return passed


# ──────────────────────────────────────────────────────────
# SDF プロパティ用ヘルパー
# ──────────────────────────────────────────────────────────

def get_synth_props_for_sdf(candidate: dict) -> dict:
    """SDF 保存時に SetProp する用のキー/値ペアを返す"""
    return {
        "SA_Score":      str(candidate.get("SA_Score", "")),
        "QED":           str(candidate.get("QED", "")),
        "PAINS_OK":      str(candidate.get("PAINS_OK", "")),
        "PAINS_Alerts":  candidate.get("PAINS_Alerts", ""),
        "BRENK_OK":      str(candidate.get("BRENK_OK", "")),
        "BRENK_Alerts":  candidate.get("BRENK_Alerts", ""),
        "Synth_Pass":    str(candidate.get("Synth_Pass", "")),
    }


# ──────────────────────────────────────────────────────────
# SMILES から全スコアを計算するユーティリティ
# ──────────────────────────────────────────────────────────

def score_from_smiles(smiles: str) -> dict:
    """SMILES 文字列から全合成容易性スコアを計算して dict で返す"""
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    sa = calc_sa_score(mol)
    qed = calc_qed(mol)
    pains_ok, pains_alerts = check_pains(mol)
    brenk_ok, brenk_alerts = check_brenk(mol)
    sa_ok = (not math.isnan(sa)) and sa <= 6.0
    return {
        "sa_score":      sa,
        "qed":           qed,
        "pains_ok":      pains_ok,
        "brenk_ok":      brenk_ok,
        "pains_alerts":  "; ".join(pains_alerts),
        "brenk_alerts":  "; ".join(brenk_alerts),
        "synth_pass":    sa_ok and pains_ok and brenk_ok,
    }


# ──────────────────────────────────────────────────────────
# 動作テスト
# ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    test_smiles = [
        ("aspirin",    "CC(=O)Oc1ccccc1C(=O)O"),
        ("caffeine",   "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
        ("benzene",    "c1ccccc1"),
        ("ibuprofen",  "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("taxol_frag", "CC1(C)CC(O)C(=O)C(C)=C1"),
    ]

    print("=" * 70)
    print(f"  {'名前':<15} {'SA':>5} {'QED':>6} {'PAINS':>6} {'BRENK':>6} {'合格':>4}")
    print("=" * 70)

    for name, smi in test_smiles:
        mol = Chem.MolFromSmiles(smi)
        sa  = calc_sa_score(mol)
        qed = calc_qed(mol)
        pains_ok, _ = check_pains(mol)
        brenk_ok, _ = check_brenk(mol)
        sa_ok = (not math.isnan(sa)) and sa <= 6.0
        synth = sa_ok and pains_ok and brenk_ok
        print(f"  {name:<15} {sa:>5.2f} {qed:>6.3f} "
              f"{'OK' if pains_ok else 'NG':>6} "
              f"{'OK' if brenk_ok else 'NG':>6} "
              f"{'OK' if synth else 'NG':>4}")

    print("=" * 70)
