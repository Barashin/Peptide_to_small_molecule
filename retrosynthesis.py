#!/usr/bin/env python3
"""
retrosynthesis.py
=================
候補分子の逆合成解析モジュール。

3 段階のアプローチ:

  1. BRICS 逆合成 (RDKit)
     - BRICS ルール (合成的に意味のある結合切断) でフラグメントに分解
     - 各フラグメントの SA Score を計算 → 市販品レベルか判定
     - 切断箇所数 = 合成ステップ数の目安

  2. RECAP 逆合成 (RDKit)
     - RECAP (Retrosynthetic Combinatorial Analysis Procedure)
     - メディシナルケミストリーで一般的な反応で切断
     - BRICS と相補的: 異なるルールセットで冗長性を確保

  3. ASKCOS API (オプション, MIT)
     - 外部 API を呼んで AI ベースの逆合成ルートを取得
     - URL が設定されている場合のみ使用
"""

import json
import math
import os
from pathlib import Path
from typing import Optional

from rdkit import Chem
from rdkit.Chem import BRICS, Recap, Descriptors, AllChem


# ──────────────────────────────────────────────────────────
# SA Score (reuse from synthesizability.py)
# ──────────────────────────────────────────────────────────

def _sa_score(mol) -> float:
    try:
        from synthesizability import calc_sa_score
        return calc_sa_score(mol)
    except ImportError:
        return float("nan")


# ──────────────────────────────────────────────────────────
# BRICS 反応タイプの説明 (日本語)
# ──────────────────────────────────────────────────────────

BRICS_BOND_TYPES = {
    (1,):  "C-C 単結合 (sp3)",
    (3,):  "エーテル / エステル C-O",
    (5,):  "アミド C-N",
    (6,):  "芳香族 C=C (ビアリール等)",
    (7,):  "C-N (アミン / スルホンアミド)",
    (8,):  "C-O (アルキル エーテル)",
    (9,):  "C-N (ラクタム / 環状アミド)",
    (10,): "N-N (ヒドラジン / ヒドラジド)",
    (11,): "C-S (チオエーテル)",
    (12,): "C=C (オレフィン メタセシス)",
    (13,): "C-C (環外結合)",
    (14,): "芳香環 C-ヘテロ原子",
    (15,): "C-N (複素環)",
    (16,): "芳香環 C-C (Suzuki 等)",
}


def _brics_label_to_reaction(label: str) -> str:
    """BRICS ダミーアトムラベル [N*] → 反応タイプ名を返す"""
    import re
    nums = re.findall(r"\[(\d+)\*\]", label)
    if not nums:
        return "不明"
    key = tuple(sorted(set(int(n) for n in nums)))
    return BRICS_BOND_TYPES.get(key, f"BRICS type {key}")


# ──────────────────────────────────────────────────────────
# BRICS 逆合成
# ──────────────────────────────────────────────────────────

def brics_retrosynthesis(mol) -> dict:
    """
    BRICS ルールに基づく逆合成解析。

    Returns:
        {
          "n_cuts":      int,       # 切断箇所数 (≈ 合成ステップ数)
          "fragments":   [ {smiles, sa_score, n_atoms, buyable} ],
          "reactions":   [ str ],   # 推定される反応タイプ
          "max_frag_sa": float,     # 最も合成困難なフラグメントの SA
          "feasibility": str,       # "容易" / "中程度" / "困難"
        }
    """
    if mol is None:
        return {"n_cuts": 0, "fragments": [], "reactions": [],
                "max_frag_sa": float("nan"), "feasibility": "解析不可"}

    frag_smiles_set = BRICS.BRICSDecompose(mol)

    fragments = []
    reactions = set()
    for smi in frag_smiles_set:
        # 反応タイプを抽出
        reaction_type = _brics_label_to_reaction(smi)
        if reaction_type != "不明":
            reactions.add(reaction_type)

        # ダミーアトム [N*] を除去して実フラグメントに変換
        import re
        clean_smi = re.sub(r"\[\d+\*\]", "[H]", smi)

        clean_mol = Chem.MolFromSmiles(clean_smi)
        if clean_mol is None:
            continue

        n_atoms = clean_mol.GetNumHeavyAtoms()
        sa = _sa_score(clean_mol)
        # 市販品判定: 小さいフラグメントは基本的に入手可能
        buyable = (n_atoms <= 3) or (
            (not math.isnan(sa)) and sa <= 3.5 and n_atoms <= 15
        )

        fragments.append({
            "smiles":   Chem.MolToSmiles(clean_mol),
            "sa_score":  sa,
            "n_atoms":  n_atoms,
            "buyable":  buyable,
        })

    n_cuts = max(len(frag_smiles_set) - 1, 0)

    # 最大 SA Score
    sa_values = [f["sa_score"] for f in fragments if not math.isnan(f["sa_score"])]
    max_sa = max(sa_values) if sa_values else float("nan")

    # 元分子自体の SA Score
    mol_sa = _sa_score(mol)
    mol_buyable = (not math.isnan(mol_sa)) and mol_sa <= 3.0

    # 合成可能性判定
    all_buyable = all(f["buyable"] for f in fragments) if fragments else False
    if mol_buyable:
        feasibility = "市販品レベル (SA≤3.0, そのまま購入可能な可能性大)"
    elif n_cuts <= 2 and all_buyable:
        feasibility = "容易 (2 ステップ以内, 全フラグメント市販品レベル)"
    elif n_cuts <= 4 and all_buyable:
        feasibility = "比較的容易 (市販フラグメントの組み合わせ)"
    elif n_cuts <= 4 and max_sa is not None and max_sa <= 4.0:
        feasibility = "中程度 (3-4 ステップ, フラグメントは合成容易)"
    elif n_cuts <= 4:
        feasibility = "やや困難 (特殊フラグメントを含む)"
    else:
        feasibility = "困難 (5 ステップ以上の多段階合成)"

    return {
        "n_cuts":      n_cuts,
        "fragments":   fragments,
        "reactions":   sorted(reactions),
        "max_frag_sa": round(max_sa, 2) if not math.isnan(max_sa) else None,
        "feasibility": feasibility,
    }


# ──────────────────────────────────────────────────────────
# RECAP 逆合成
# ──────────────────────────────────────────────────────────

RECAP_REACTION_NAMES = {
    "amide":         "アミド結合形成",
    "ester":         "エステル化",
    "amine":         "還元的アミノ化",
    "urea":          "ウレア形成",
    "ether":         "ウィリアムソンエーテル合成",
    "olefin":        "オレフィンメタセシス",
    "aromatic_n":    "芳香族 N-アルキル化",
    "lactam_n":      "ラクタム N-アルキル化",
    "aromatic_c":    "Suzuki / Heck カップリング",
    "sulfonamide":   "スルホンアミド形成",
    "thioether":     "チオエーテル形成",
}


def recap_retrosynthesis(mol) -> dict:
    """
    RECAP ルールに基づく逆合成解析。
    BRICS より保守的な切断ルールセット。

    Returns:
        {
          "n_children":  int,       # 直接の子ノード数
          "all_leaves":  [ {smiles, sa_score, n_atoms} ],
          "tree_depth":  int,       # 再帰的分解の最大深さ
        }
    """
    if mol is None:
        return {"n_children": 0, "all_leaves": [], "tree_depth": 0}

    tree = Recap.RecapDecompose(mol)

    # 再帰的に全リーフノードを収集
    def _get_leaves(node, depth=0):
        if not node.children:
            return [(node.smiles, depth)]
        results = []
        for child_smi, child_node in node.children.items():
            results.extend(_get_leaves(child_node, depth + 1))
        return results

    leaf_data = _get_leaves(tree)
    max_depth = max(d for _, d in leaf_data) if leaf_data else 0

    leaves = []
    seen = set()
    for smi, _ in leaf_data:
        if smi in seen:
            continue
        seen.add(smi)
        # RECAP は [*] をダミーアトムとして使う
        clean_smi = smi.replace("*", "[H]")
        clean_mol = Chem.MolFromSmiles(clean_smi)
        if clean_mol is None:
            clean_mol = Chem.MolFromSmiles(smi)
        if clean_mol is None:
            continue
        sa = _sa_score(clean_mol)
        leaves.append({
            "smiles":  Chem.MolToSmiles(clean_mol),
            "sa_score": sa,
            "n_atoms": clean_mol.GetNumHeavyAtoms(),
        })

    return {
        "n_children":  len(tree.children),
        "all_leaves":  leaves,
        "tree_depth":  max_depth,
    }


# ──────────────────────────────────────────────────────────
# ASKCOS API (オプション)
# ──────────────────────────────────────────────────────────

ASKCOS_DEFAULT_URL = "https://askcos.mit.edu/api/v2"


def askcos_retrosynthesis(
    smiles: str,
    api_url: str = ASKCOS_DEFAULT_URL,
    max_depth: int = 5,
    timeout: int = 30,
) -> Optional[dict]:
    """
    ASKCOS (MIT) API を使った AI ベースの逆合成解析。
    API が利用不可の場合は None を返す。

    Returns:
        {
          "status":     str,
          "n_routes":   int,
          "best_route": {
              "n_steps":    int,
              "reactions":  [ {smiles, template} ],
              "precursors": [ {smiles, buyable} ],
          },
          "all_routes": [ ... ],
        }
    """
    try:
        import urllib.request
        import urllib.error

        payload = json.dumps({
            "smiles": smiles,
            "max_depth": max_depth,
            "max_branching": 25,
            "expansion_time": 60,
        }).encode("utf-8")

        req = urllib.request.Request(
            f"{api_url}/tree-search/",
            data=payload,
            headers={"Content-Type": "application/json"},
        )

        with urllib.request.urlopen(req, timeout=timeout) as resp:
            result = json.loads(resp.read().decode("utf-8"))

        if not result.get("result"):
            return {"status": "ルートなし", "n_routes": 0,
                    "best_route": None, "all_routes": []}

        routes = result["result"]
        best = routes[0] if routes else None
        return {
            "status":     "成功",
            "n_routes":   len(routes),
            "best_route": best,
            "all_routes": routes,
        }

    except (urllib.error.URLError, urllib.error.HTTPError, OSError) as e:
        return None
    except Exception:
        return None


# ──────────────────────────────────────────────────────────
# AiZynthFinder (オプション, ローカル実行)
# ──────────────────────────────────────────────────────────

_AIZF_FINDER = None  # シングルトン (初回のみ初期化)


def _get_aizf_config() -> str:
    """AiZynthFinder config.yml のパスを自動検索"""
    candidates = [
        Path(__file__).resolve().parent / "aizynthfinder_data" / "config.yml",
        Path(__file__).resolve().parent.parent / "aizynthfinder_data" / "config.yml",
        Path.home() / "aizynthfinder_data" / "config.yml",
    ]
    for p in candidates:
        if p.is_file():
            return str(p)
    return ""


def _extract_route(node: dict) -> dict:
    """
    AiZynthFinder のルート dict (木構造) から
    反応・中間体・出発物を再帰的に抽出する。

    Returns:
        {
          "reactions": [
            {
              "smarts": str,
              "product": str,
              "reactants": [str, ...],
              "classification": str,
              "template": str,
            }
          ],
          "starting_materials": [
            {"smiles": str, "in_stock": bool}
          ],
        }
    """
    reactions = []
    starting_materials = []

    def _traverse(n):
        if n.get("type") == "mol" or n.get("is_chemical"):
            if n.get("in_stock") and not n.get("children"):
                starting_materials.append({
                    "smiles": n["smiles"],
                    "in_stock": True,
                })
            elif n.get("children"):
                for child in n["children"]:
                    _traverse(child)
            else:
                # 在庫なし & 子なし (leaf but not in stock)
                starting_materials.append({
                    "smiles": n["smiles"],
                    "in_stock": False,
                })
        elif n.get("is_reaction") or n.get("type") == "reaction":
            reactant_smiles = []
            parent_smiles = ""
            # 親 mol = reaction の predecessor (product)
            # AiZynthFinder dict では reaction の children = reactants
            if n.get("children"):
                for child in n["children"]:
                    if child.get("is_chemical") or child.get("type") == "mol":
                        reactant_smiles.append(child["smiles"])
                    _traverse(child)

            meta = n.get("metadata", {})
            reactions.append({
                "smarts": n.get("smiles", ""),
                "reactants": reactant_smiles,
                "classification": meta.get("classification", ""),
                "template": meta.get("template", ""),
            })

    _traverse(node)
    return {
        "reactions": reactions,
        "starting_materials": starting_materials,
    }


def _route_to_forward_steps(route_dict: dict) -> list[dict]:
    """
    AiZynthFinder の逆合成ルート (ターゲット→原料) を
    前向き合成シーケンス (原料→ターゲット) に変換。

    Returns:
        [
          {"step": 1,
           "substrate_smi": "...",
           "reagent_smi": "..." or None,
           "product_smi": "...",
           "reaction_en": "...",
           "reagents": "...",
           "solvent": "...",
           "temp": "...",
           "yield_range": "...",
           "classification": "...",
          },
          ...
        ]
    """
    # 木構造を逆順に辿って前向きステップを生成
    steps = []

    def _traverse(node) -> str | None:
        """ノードを辿り、このノードの SMILES を返す。
        反応ノードに遭遇したら step を記録する。"""
        if node.get("type") == "mol" or node.get("is_chemical"):
            if not node.get("children"):
                # leaf: starting material or unsolved
                return node.get("smiles")
            # mol with children → has a reaction
            for child in node.get("children", []):
                _traverse(child)
            return node.get("smiles")
        elif node.get("is_reaction") or node.get("type") == "reaction":
            # Gather reactant SMILES (after recursing into them)
            reactant_smiles = []
            for child in node.get("children", []):
                smi = _traverse(child)
                if smi:
                    reactant_smiles.append(smi)

            # This reaction's product is the parent mol node
            # We'll set it below after we know all products
            meta = node.get("metadata", {})
            cls = meta.get("classification", "")

            steps.append({
                "reactants": reactant_smiles,
                "classification": cls,
                "template": meta.get("template", ""),
                "smarts": node.get("smiles", ""),
            })
            return None
        return None

    _traverse(route_dict)

    # 木構造からproductを決定: 各反応の生成物を特定する
    # AiZynthFinder の木: mol → reaction → [reactant_mol, ...]
    # 前向き: reactants → product (= parent mol)
    def _assign_products(node):
        """各反応の生成物を親molのSMILESで埋める"""
        if (node.get("type") == "mol" or node.get("is_chemical")) and node.get("children"):
            for child in node["children"]:
                if child.get("is_reaction") or child.get("type") == "reaction":
                    child["_product_smi"] = node["smiles"]
                    for grandchild in child.get("children", []):
                        _assign_products(grandchild)

    _assign_products(route_dict)

    # 再度走査して product を取得
    steps2 = []

    def _traverse2(node):
        if (node.get("type") == "mol" or node.get("is_chemical")) and node.get("children"):
            for child in node["children"]:
                _traverse2(child)
        elif node.get("is_reaction") or node.get("type") == "reaction":
            for child in node.get("children", []):
                _traverse2(child)
            meta = node.get("metadata", {})
            reactant_smiles = [c["smiles"] for c in node.get("children", [])
                               if c.get("type") == "mol" or c.get("is_chemical")]
            steps2.append({
                "reactants": reactant_smiles,
                "product_smi": node.get("_product_smi", ""),
                "classification": meta.get("classification", ""),
                "template": meta.get("template", ""),
                "smarts": node.get("smiles", ""),
            })

    _traverse2(route_dict)

    # 前向き順序: deepest reactions first (= already in order from DFS post-order)
    forward_steps = []
    for i, s in enumerate(steps2):
        reactants = s["reactants"]
        # substrate = first (or largest) reactant, reagent = rest
        substrate_smi = reactants[0] if reactants else ""
        reagent_smi = reactants[1] if len(reactants) > 1 else None

        forward_steps.append({
            "step": i + 1,
            "substrate_smi": substrate_smi,
            "reagent_smi": reagent_smi,
            "product_smi": s["product_smi"],
            "reaction_en": _classify_reaction(s["classification"], s["smarts"]),
            "reagents": "",
            "solvent": "",
            "temp": "",
            "yield_range": "",
            "classification": s["classification"],
        })

    return forward_steps


def _classify_reaction(classification: str, smarts: str) -> str:
    """反応分類から英語の反応名を推定する"""
    cls_lower = classification.lower() if classification else ""

    # NameRXN classification patterns
    if "negishi" in cls_lower:
        return "Negishi Coupling"
    if "suzuki" in cls_lower:
        return "Suzuki Coupling"
    if "heck" in cls_lower:
        return "Heck Reaction"
    if "buchwald" in cls_lower or "hartwig" in cls_lower:
        return "Buchwald-Hartwig Amination"
    if "sonogashira" in cls_lower:
        return "Sonogashira Coupling"
    if "grignard" in cls_lower:
        return "Grignard Reaction"
    if "wittig" in cls_lower:
        return "Wittig Reaction"
    if "aldol" in cls_lower:
        return "Aldol Reaction"
    if "amide" in cls_lower:
        return "Amide Coupling"
    if "reduction" in cls_lower:
        return "Reduction"
    if "oxidation" in cls_lower:
        return "Oxidation"
    if "hydrogenation" in cls_lower:
        return "Hydrogenation"
    if "deprotection" in cls_lower:
        return "Deprotection"
    if "protection" in cls_lower:
        return "Protection"
    if "alkylation" in cls_lower:
        return "Alkylation"
    if "acylation" in cls_lower:
        return "Acylation"
    if "ester" in cls_lower:
        return "Esterification"
    if "cyclization" in cls_lower:
        return "Cyclization"

    # SMARTS pattern analysis for "0.0 Unrecognized" cases
    if smarts:
        if ">>" in smarts:
            prod_side, react_side = smarts.split(">>", 1)
            if "O=" in react_side and "O=" not in prod_side:
                return "Reduction (C=O to CH2)"
            if "[OH]" in react_side or "O=" in react_side:
                if "N" in prod_side:
                    return "Amide Coupling"
            if "Br" in react_side or "Cl" in react_side or "I" in react_side:
                if "[Zn" in react_side:
                    return "Negishi Coupling"
                if "B(" in react_side or "B]" in react_side:
                    return "Suzuki Coupling"
                return "Cross-Coupling"
            if "N" in react_side and "C(=O)" in prod_side:
                return "Amide Coupling"

    if classification and "unrecognized" not in cls_lower:
        return classification

    return "Synthetic Transformation"


def aizynthfinder_retrosynthesis(
    smiles: str,
    config_path: str = "",
    time_limit: int = 120,
    iteration_limit: int = 100,
    max_transforms: int = 6,
) -> dict | None:
    """
    AiZynthFinder で逆合成解析を実行。

    Returns:
        {
          "status": str,           # "成功" / "ルートなし" / "エラー"
          "n_routes": int,
          "routes": [
            {
              "score": float,
              "n_steps": int,
              "reactions": [...],
              "starting_materials": [...],
              "forward_steps": [...],  # 前向き合成シーケンス
            }
          ],
        }
    or None if AiZynthFinder is not installed.
    """
    global _AIZF_FINDER

    try:
        from aizynthfinder.aizynthfinder import AiZynthFinder
    except ImportError:
        return None

    # 設定ファイル
    if not config_path:
        config_path = _get_aizf_config()
    if not config_path:
        return {"status": "エラー (config.yml が見つからない)",
                "n_routes": 0, "routes": []}

    try:
        # Finder の初期化 (シングルトン)
        if _AIZF_FINDER is None:
            import warnings
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                _AIZF_FINDER = AiZynthFinder(configfile=config_path)
                _AIZF_FINDER.expansion_policy.select("uspto")
                _AIZF_FINDER.filter_policy.select("uspto")
                _AIZF_FINDER.stock.select("zinc")

        finder = _AIZF_FINDER
        finder.config.search.time_limit = time_limit
        finder.config.search.iteration_limit = iteration_limit
        finder.config.search.max_transforms = max_transforms
        finder.target_smiles = smiles

        finder.tree_search()
        finder.build_routes()

        rc = finder.routes
        if len(rc) == 0:
            return {"status": "ルートなし", "n_routes": 0, "routes": []}

        rc.make_dicts()
        scores = rc.scores  # [{"state score": float}, ...]

        routes_out = []
        for i, route_dict in enumerate(rc.dicts):
            score_val = list(scores[i].values())[0] if i < len(scores) else 0.0

            extracted = _extract_route(route_dict)
            forward_steps = _route_to_forward_steps(route_dict)

            routes_out.append({
                "score": round(score_val, 4),
                "n_steps": len(extracted["reactions"]),
                "reactions": extracted["reactions"],
                "starting_materials": extracted["starting_materials"],
                "forward_steps": forward_steps,
                "raw_dict": route_dict,
            })

        # スコア降順でソート
        routes_out.sort(key=lambda r: -r["score"])

        return {
            "status": "成功",
            "n_routes": len(routes_out),
            "routes": routes_out,
        }

    except Exception as e:
        return {"status": f"エラー ({e})", "n_routes": 0, "routes": []}


# ──────────────────────────────────────────────────────────
# 統合解析関数
# ──────────────────────────────────────────────────────────

def analyze_retrosynthesis(
    smiles: str,
    name: str = "",
    use_askcos: bool = False,
    askcos_url: str = ASKCOS_DEFAULT_URL,
    use_aizynthfinder: bool = False,
    aizynthfinder_config: str = "",
    verbose: bool = True,
) -> dict:
    """
    1 分子に対する統合的な逆合成解析。

    Returns:
        {
          "name":     str,
          "smiles":   str,
          "brics":    dict,   # BRICS 解析結果
          "recap":    dict,   # RECAP 解析結果
          "askcos":   dict | None,  # ASKCOS 結果 (API 使用時のみ)
          "aizynthfinder": dict | None,  # AiZynthFinder 結果
          "summary":  str,    # 人間可読なサマリー
          "overall_feasibility": str,  # 総合判定
        }
    """
    mol = Chem.MolFromSmiles(smiles) if smiles else None

    brics_result = brics_retrosynthesis(mol)
    recap_result = recap_retrosynthesis(mol)

    askcos_result = None
    if use_askcos:
        askcos_result = askcos_retrosynthesis(smiles, api_url=askcos_url)

    aizf_result = None
    if use_aizynthfinder:
        aizf_result = aizynthfinder_retrosynthesis(
            smiles, config_path=aizynthfinder_config)

    # 総合判定
    overall = _overall_feasibility(
        brics_result, recap_result, askcos_result, aizf_result)

    # サマリー生成
    summary = _format_summary(name, smiles, brics_result, recap_result,
                               askcos_result, overall, aizf_result)
    if verbose:
        print(summary)

    return {
        "name":     name,
        "smiles":   smiles,
        "brics":    brics_result,
        "recap":    recap_result,
        "askcos":   askcos_result,
        "aizynthfinder": aizf_result,
        "summary":  summary,
        "overall_feasibility": overall,
    }


def _overall_feasibility(
    brics: dict, recap: dict,
    askcos: Optional[dict],
    aizf: Optional[dict] = None,
) -> str:
    """BRICS/RECAP の結果から総合的な合成可能性を判定"""
    n_cuts = brics.get("n_cuts", 0)
    max_sa = brics.get("max_frag_sa")
    all_buyable = all(f.get("buyable", False) for f in brics.get("fragments", []))
    recap_depth = recap.get("tree_depth", 0)

    # AiZynthFinder が使えてルートが見つかった場合
    if aizf and aizf.get("n_routes", 0) > 0:
        best = aizf["routes"][0]
        all_in_stock = all(
            sm.get("in_stock") for sm in best.get("starting_materials", []))
        n_steps = best.get("n_steps", 0)
        if all_in_stock and n_steps <= 2:
            return "容易 (AI ルートあり, 全原料市販品)"
        if all_in_stock:
            return f"合成可能 (AI {n_steps} ステップ, 全原料市販品)"
        return f"合成可能 (AI {n_steps} ステップ)"

    # ASKCOS が使えてルートが見つかった場合
    if askcos and askcos.get("n_routes", 0) > 0:
        return "合成可能 (AI 逆合成ルートあり)"

    # BRICS + RECAP の組み合わせ判定
    if n_cuts == 0:
        return "市販品レベル (切断不要)"
    if all_buyable and n_cuts <= 2:
        return "容易 (1-2 ステップ, 市販フラグメントの組み合わせ)"
    if all_buyable and n_cuts <= 4:
        return "比較的容易 (市販フラグメントの組み合わせ)"
    if n_cuts <= 3 and max_sa is not None and max_sa <= 4.0:
        return "中程度 (標準的な有機合成で到達可能)"
    if n_cuts <= 4 and max_sa is not None and max_sa <= 5.0:
        return "中程度 (やや複雑だが到達可能)"
    if n_cuts <= 5:
        return "やや困難 (特殊試薬 or 多段階合成)"
    return "困難 (6 ステップ以上の多段階合成, 要専門家相談)"


def _format_summary(
    name: str, smiles: str,
    brics: dict, recap: dict,
    askcos: Optional[dict], overall: str,
    aizf: Optional[dict] = None,
) -> str:
    """人間可読な逆合成サマリーを生成"""
    lines = []
    label = name or smiles[:40]
    lines.append(f"\n  [{label}] 逆合成解析")
    lines.append(f"  {'─' * 50}")

    # BRICS
    lines.append(f"  BRICS 切断: {brics['n_cuts']} 箇所")
    if brics["reactions"]:
        lines.append(f"    推定反応: {', '.join(brics['reactions'])}")
    for i, f in enumerate(brics.get("fragments", []), 1):
        sa_s = f"{f['sa_score']:.1f}" if not math.isnan(f["sa_score"]) else "?"
        buy = "市販品" if f["buyable"] else "要合成"
        lines.append(f"    フラグメント {i}: {f['smiles']:<30} "
                     f"(SA={sa_s}, {f['n_atoms']}原子, {buy})")

    # RECAP
    lines.append(f"  RECAP 分解深度: {recap['tree_depth']}, "
                 f"リーフ数: {len(recap['all_leaves'])}")

    # ASKCOS
    if askcos is not None:
        if askcos.get("n_routes", 0) > 0:
            lines.append(f"  ASKCOS: {askcos['n_routes']} ルート発見")
        else:
            lines.append(f"  ASKCOS: {askcos.get('status', 'N/A')}")

    # AiZynthFinder
    if aizf is not None:
        if aizf.get("n_routes", 0) > 0:
            best = aizf["routes"][0]
            n_sm = len(best.get("starting_materials", []))
            n_in = sum(1 for s in best.get("starting_materials", [])
                       if s.get("in_stock"))
            lines.append(f"  AiZynthFinder: {aizf['n_routes']} ルート発見 "
                         f"(best: {best['n_steps']} steps, "
                         f"score={best['score']:.3f}, "
                         f"原料 {n_in}/{n_sm} 市販品)")
        else:
            lines.append(f"  AiZynthFinder: {aizf.get('status', 'N/A')}")

    lines.append(f"  総合判定: {overall}")
    return "\n".join(lines)


# ──────────────────────────────────────────────────────────
# バッチ解析
# ──────────────────────────────────────────────────────────

def batch_retrosynthesis(
    candidates: list,
    use_askcos: bool = False,
    askcos_url: str = ASKCOS_DEFAULT_URL,
    use_aizynthfinder: bool = False,
    aizynthfinder_config: str = "",
    verbose: bool = True,
) -> list:
    """
    候補分子リストに対するバッチ逆合成解析。

    Args:
        candidates: [{"name": str, "smiles": str, ...}, ...]

    Returns:
        [{"name", "smiles", "brics", "recap", "overall_feasibility", ...}, ...]
    """
    results = []
    for c in candidates:
        smiles = c.get("smiles", "")
        name = c.get("name", "")
        if not smiles:
            continue
        r = analyze_retrosynthesis(
            smiles, name=name,
            use_askcos=use_askcos, askcos_url=askcos_url,
            use_aizynthfinder=use_aizynthfinder,
            aizynthfinder_config=aizynthfinder_config,
            verbose=verbose,
        )
        results.append(r)

    if verbose and results:
        print(f"\n  {'=' * 60}")
        print(f"  逆合成解析サマリー: {len(results)} 分子")
        print(f"  {'=' * 60}")
        easy = sum(1 for r in results
                   if "容易" in r["overall_feasibility"]
                   or "市販品" in r["overall_feasibility"])
        medium = sum(1 for r in results
                     if "中程度" in r["overall_feasibility"])
        hard = sum(1 for r in results
                   if "困難" in r["overall_feasibility"])
        print(f"    容易:   {easy}/{len(results)}")
        print(f"    中程度: {medium}/{len(results)}")
        print(f"    困難:   {hard}/{len(results)}")

    return results


# ──────────────────────────────────────────────────────────
# レポート出力
# ──────────────────────────────────────────────────────────

def save_retrosynthesis_report(
    results: list,
    output_dir: str = "results",
    filename: str = "retrosynthesis_report.json",
):
    """逆合成解析結果を JSON で保存"""
    out_path = Path(output_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # JSON シリアライズ可能にする
    serializable = []
    for r in results:
        entry = {
            "name":   r["name"],
            "smiles": r["smiles"],
            "overall_feasibility": r["overall_feasibility"],
            "brics": {
                "n_cuts":      r["brics"]["n_cuts"],
                "max_frag_sa": r["brics"]["max_frag_sa"],
                "feasibility": r["brics"]["feasibility"],
                "fragments":   r["brics"]["fragments"],
                "reactions":   r["brics"]["reactions"],
            },
            "recap": {
                "n_children": r["recap"]["n_children"],
                "tree_depth": r["recap"]["tree_depth"],
                "all_leaves": r["recap"]["all_leaves"],
            },
        }
        if r.get("askcos"):
            entry["askcos"] = {
                "status":   r["askcos"].get("status"),
                "n_routes": r["askcos"].get("n_routes"),
            }
        if r.get("aizynthfinder"):
            aizf = r["aizynthfinder"]
            entry["aizynthfinder"] = {
                "status":   aizf.get("status"),
                "n_routes": aizf.get("n_routes"),
            }
            if aizf.get("routes"):
                best = aizf["routes"][0]
                entry["aizynthfinder"]["best_route"] = {
                    "score": best.get("score"),
                    "n_steps": best.get("n_steps"),
                    "starting_materials": best.get("starting_materials"),
                }
        serializable.append(entry)

    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(serializable, f, indent=2, ensure_ascii=False)

    print(f"\n  逆合成レポート保存: {out_path}")
    return str(out_path)


def save_retrosynthesis_csv(
    results: list,
    output_dir: str = "results",
    filename: str = "retrosynthesis_summary.csv",
):
    """逆合成解析結果を CSV で保存"""
    import csv
    out_path = Path(output_dir) / filename
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fields = [
        "name", "smiles", "overall_feasibility",
        "brics_n_cuts", "brics_max_frag_sa", "brics_feasibility",
        "brics_reactions", "brics_n_buyable", "brics_n_fragments",
        "recap_tree_depth", "recap_n_leaves",
    ]

    with open(out_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in results:
            frags = r["brics"].get("fragments", [])
            w.writerow({
                "name":                  r["name"],
                "smiles":                r["smiles"],
                "overall_feasibility":   r["overall_feasibility"],
                "brics_n_cuts":          r["brics"]["n_cuts"],
                "brics_max_frag_sa":     r["brics"]["max_frag_sa"],
                "brics_feasibility":     r["brics"]["feasibility"],
                "brics_reactions":       "; ".join(r["brics"]["reactions"]),
                "brics_n_buyable":       sum(1 for f in frags if f.get("buyable")),
                "brics_n_fragments":     len(frags),
                "recap_tree_depth":      r["recap"]["tree_depth"],
                "recap_n_leaves":        len(r["recap"]["all_leaves"]),
            })

    print(f"  逆合成 CSV 保存: {out_path}")
    return str(out_path)


# ──────────────────────────────────────────────────────────
# テスト
# ──────────────────────────────────────────────────────────

if __name__ == "__main__":
    test_molecules = [
        ("aspirin",    "CC(=O)Oc1ccccc1C(=O)O"),
        ("ibuprofen",  "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
        ("caffeine",   "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
        ("indole_deriv", "CCc1cccc2[nH]ccc12"),
    ]

    all_results = []
    for name, smi in test_molecules:
        r = analyze_retrosynthesis(smi, name=name, verbose=True)
        all_results.append(r)

    save_retrosynthesis_report(all_results, output_dir=".")
    save_retrosynthesis_csv(all_results, output_dir=".")
