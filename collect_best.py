#!/usr/bin/env python3
"""
collect_best.py
===============
全ドッキング結果を集約し、結合力・Ligand Efficiency の上位分子を
Result_Best/ ディレクトリにまとめる。

選抜基準:
  1. dock_score あり + ドッキング SDF が実在
  2. DrugLike = True (bridge 結果のみ判定可; peptidomimetics は常に含める)
  3. smina スコア上位 + LE 上位の 2 軸で統合ランキング → Top 5 選抜

出力:
  Result_Best/
  ├── <rank>_<score>_<name>.sdf        # ドッキングポーズ SDF
  ├── retrosynthesis/                  # 逆合成解析
  │   ├── 01_retrosynthesis_*.png      # 合成スキーム図
  │   └── retrosynthesis_detail.json   # 詳細データ
  ├── summary.csv                      # 全選抜分子の詳細
  └── README.md                        # 見方と考察
"""

import os
import csv
import shutil
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors, AllChem

try:
    from synthesizability import score_from_smiles as _score_from_smiles
    _HAS_SYNTH = True
except ImportError:
    _HAS_SYNTH = False

try:
    from retrosynthesis import analyze_retrosynthesis as _retro_analyze
    _HAS_RETRO = True
except ImportError:
    _HAS_RETRO = False

try:
    from retrosynthesis import aizynthfinder_retrosynthesis as _aizf_retro
    _HAS_AIZF = True
except ImportError:
    _HAS_AIZF = False

try:
    from rdkit.Chem import Draw
    from PIL import Image, ImageDraw, ImageFont
    _HAS_DRAW = True
except ImportError:
    _HAS_DRAW = False


# ──────────────────────────────────────────────────────────
# 設定
# ──────────────────────────────────────────────────────────

RESULTS_DIR  = Path("results")
OUT_DIR      = Path("Result_Best")
TOP_N        = 5    # 最終選抜数 (スコア + LE 統合 Top N)
SCORE_CUTOFF = -5.0 # これより弱い結合は除外


# ──────────────────────────────────────────────────────────
# ユーティリティ
# ──────────────────────────────────────────────────────────

from utils.ligand_efficiency import calc_hac, calc_le, le_grade


def _enrich_synth(record: dict) -> dict:
    """SMILES から合成容易性スコアを計算して record に追加"""
    if not _HAS_SYNTH:
        return record
    smiles = record.get("smiles", "")
    if not smiles:
        return record
    scores = _score_from_smiles(smiles)
    record["sa_score"]  = scores["sa_score"]
    record["qed"]       = scores["qed"]
    record["pains_ok"]  = scores["pains_ok"]
    record["brenk_ok"]  = scores["brenk_ok"]
    return record


def _enrich_retro(record: dict) -> dict:
    """SMILES から逆合成解析を追加"""
    if not _HAS_RETRO:
        return record
    smiles = record.get("smiles", "")
    if not smiles:
        return record
    r = _retro_analyze(smiles, name=record.get("name", ""), verbose=False)
    record["retro_feasibility"] = r["overall_feasibility"]
    record["retro_n_cuts"]      = r["brics"]["n_cuts"]
    record["retro_reactions"]   = "; ".join(r["brics"]["reactions"])
    return record


# ──────────────────────────────────────────────────────────
# 結果読み込み
# ──────────────────────────────────────────────────────────

def load_all_results() -> list[dict]:
    records = []

    # ── ペプチドミメティクス (STEP 3) ──
    csv_path = RESULTS_DIR / "docking" / "docking_results.csv"
    if csv_path.exists():
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                score_str = row.get("best_score_kcal/mol", "N/A")
                if score_str == "N/A":
                    continue
                try:
                    score   = float(score_str)
                    smiles  = row.get("smiles", "")
                    sdf_path = row.get("output_sdf", "")
                    hac     = calc_hac(smiles)
                    rec = {
                        "name"      : row["name"],
                        "score"     : score,
                        "smiles"    : smiles,
                        "hac"       : hac,
                        "le"        : calc_le(score, hac),
                        "sdf_path"  : sdf_path,
                        "druglike"  : True,   # peptidomimetics は全件 DrugLike 済みとみなす
                        "source"    : "ペプチドミメティクス (STEP3)",
                        "pair"      : "TRP6-ALA7-ASP4",
                        "mw"        : "",
                        "logp"      : "",
                    }
                    _enrich_synth(rec)
                    _enrich_retro(rec)
                    records.append(rec)
                except ValueError:
                    pass

    # ── ファーマコフォアブリッジ (STEP 7) ──
    bridge_dir = RESULTS_DIR / "bridge"
    if bridge_dir.is_dir():
        for csv_path in sorted(bridge_dir.glob("bridge_*_results.csv")):
            with open(csv_path) as f:
                for row in csv.DictReader(f):
                    score_str = row.get("dock_score", "").strip()
                    if not score_str:
                        continue
                    try:
                        score    = float(score_str)
                        smiles   = row.get("smiles", "")
                        sdf_path = row.get("sdf_docked", "").strip()
                        druglike = row.get("DrugLike", "True").strip() == "True"
                        hac      = calc_hac(smiles)
                        rec = {
                            "name"    : row["name"],
                            "score"   : score,
                            "smiles"  : smiles,
                            "hac"     : hac,
                            "le"      : calc_le(score, hac),
                            "sdf_path": sdf_path,
                            "druglike": druglike,
                            "source"  : f"ブリッジ (STEP7) {row.get('res1','')}-{row.get('res2','')}",
                            "pair"    : f"{row.get('res1','')}-{row.get('res2','')}",
                            "mw"      : row.get("MW", ""),
                            "logp"    : row.get("LogP", ""),
                        }
                        # ブリッジCSVにSA/QEDがあればそちらを使用
                        if row.get("SA_Score"):
                            try:
                                rec["sa_score"] = float(row["SA_Score"])
                            except (ValueError, TypeError):
                                pass
                        if row.get("QED"):
                            try:
                                rec["qed"] = float(row["QED"])
                            except (ValueError, TypeError):
                                pass
                        if row.get("PAINS_OK"):
                            rec["pains_ok"] = row["PAINS_OK"].strip() == "True"
                        if row.get("BRENK_OK"):
                            rec["brenk_ok"] = row["BRENK_OK"].strip() == "True"
                        # コンフォメーション適合性 (pharmacophore_bridge も出力する)
                        if row.get("conformance_rate"):
                            try:
                                rec["conformance_rate"] = float(row["conformance_rate"])
                            except (ValueError, TypeError):
                                pass
                        if row.get("rotatable_between"):
                            try:
                                rec["rotatable_between"] = int(row["rotatable_between"])
                            except (ValueError, TypeError):
                                pass
                        # SA/QEDがまだなければSMILESから計算
                        if "sa_score" not in rec:
                            _enrich_synth(rec)
                        _enrich_retro(rec)
                        records.append(rec)
                    except ValueError:
                        pass

    # ── 剛直スキャフォールドブリッジ (STEP 7b) ──
    rigid_dir = RESULTS_DIR / "rigid_scaffold"
    if rigid_dir.is_dir():
        for csv_path in sorted(rigid_dir.glob("rigid_*_results.csv")):
            with open(csv_path) as f:
                for row in csv.DictReader(f):
                    score_str = row.get("dock_score", "").strip()
                    if not score_str:
                        continue
                    try:
                        score    = float(score_str)
                        smiles   = row.get("smiles", "")
                        sdf_path = row.get("sdf_docked", "").strip()
                        druglike = row.get("DrugLike", "True").strip() == "True"
                        hac      = calc_hac(smiles)
                        rec = {
                            "name"    : row["name"],
                            "score"   : score,
                            "smiles"  : smiles,
                            "hac"     : hac,
                            "le"      : calc_le(score, hac),
                            "sdf_path": sdf_path,
                            "druglike": druglike,
                            "source"  : f"剛直スキャフォールド (STEP7b) {row.get('res1','')}-{row.get('res2','')}",
                            "pair"    : f"{row.get('res1','')}-{row.get('res2','')}",
                            "mw"      : row.get("MW", ""),
                            "logp"    : row.get("LogP", ""),
                        }
                        # SA/QED/PAINS/BRENK
                        if row.get("SA_Score"):
                            try:
                                rec["sa_score"] = float(row["SA_Score"])
                            except (ValueError, TypeError):
                                pass
                        if row.get("QED"):
                            try:
                                rec["qed"] = float(row["QED"])
                            except (ValueError, TypeError):
                                pass
                        if row.get("PAINS_OK"):
                            rec["pains_ok"] = row["PAINS_OK"].strip() == "True"
                        if row.get("BRENK_OK"):
                            rec["brenk_ok"] = row["BRENK_OK"].strip() == "True"
                        # コンフォメーション適合性
                        if row.get("conformance_rate"):
                            try:
                                rec["conformance_rate"] = float(row["conformance_rate"])
                            except (ValueError, TypeError):
                                pass
                        if row.get("rotatable_between"):
                            try:
                                rec["rotatable_between"] = int(row["rotatable_between"])
                            except (ValueError, TypeError):
                                pass
                        if "sa_score" not in rec:
                            _enrich_synth(rec)
                        _enrich_retro(rec)
                        records.append(rec)
                    except ValueError:
                        pass

    return records


# ──────────────────────────────────────────────────────────
# 選抜ロジック
# ──────────────────────────────────────────────────────────

def select_best(records: list[dict]) -> list[dict]:
    # ① SDF が実在 + DrugLike + スコアカットオフ
    valid = [
        r for r in records
        if r["score"] <= SCORE_CUTOFF
        and r["druglike"]
        and r["sdf_path"]
        and Path(r["sdf_path"]).exists()
    ]

    if not valid:
        print("  [警告] 有効な SDF を持つ DrugLike 分子が見つかりません")
        return []

    # ② スコア順でソートし、Top N を選抜
    by_score = sorted(valid, key=lambda r: r["score"])

    # LE 上位も候補に加え、重複除去して Top N に統合
    with_le = [r for r in valid if r["le"] is not None]
    by_le   = sorted(with_le, key=lambda r: r["le"], reverse=True)

    seen  = set()
    final = []
    for r in by_score + by_le:
        if r["name"] not in seen:
            seen.add(r["name"])
            final.append(r)

    # ③ 最終スコア順でソートし Top N
    final.sort(key=lambda r: r["score"])
    return final[:TOP_N]


# ──────────────────────────────────────────────────────────
# ファイルコピーと命名
# ──────────────────────────────────────────────────────────

def safe_filename(name: str, max_len: int = 60) -> str:
    """ファイル名として安全な文字列に変換"""
    safe = name.replace(" ", "_").replace("/", "-").replace("\\", "-")
    return safe[:max_len]


def copy_sdfs(selected: list[dict], out_dir: Path) -> list[dict]:
    """選抜分子の SDF を out_dir にコピーし、コピー先パスを付与して返す"""
    enriched = []
    for rank, r in enumerate(selected, 1):
        score_tag = f"{abs(r['score']):.1f}".replace(".", "p")
        fname = f"{rank:02d}_score{score_tag}_{safe_filename(r['name'])}.sdf"
        dst   = out_dir / fname

        try:
            shutil.copy2(r["sdf_path"], dst)
            r["out_sdf"] = str(dst)
            print(f"  [{rank:02d}] {r['score']:>7.3f} kcal/mol  LE={r['le'] or 'N/A'}  → {fname}")
        except Exception as e:
            print(f"  [警告] SDF コピー失敗: {r['sdf_path']} → {e}")
            r["out_sdf"] = ""

        enriched.append(r)
    return enriched


# ──────────────────────────────────────────────────────────
# サマリー CSV
# ──────────────────────────────────────────────────────────

def write_summary_csv(selected: list[dict], out_dir: Path):
    csv_path = out_dir / "summary.csv"
    fields = ["rank", "name", "score", "hac", "le", "le_grade",
              "sa_score", "qed", "pains_ok", "brenk_ok",
              "conformance_rate", "rotatable_between",
              "retro_feasibility", "retro_n_cuts",
              "mw", "logp", "smiles", "pair", "source", "out_sdf"]
    with open(csv_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rank, r in enumerate(selected, 1):
            w.writerow({
                "rank"    : rank,
                "name"    : r["name"],
                "score"   : r["score"],
                "hac"     : r.get("hac", ""),
                "le"      : r.get("le", ""),
                "le_grade": le_grade(r.get("le")),
                "sa_score": r.get("sa_score", ""),
                "qed"     : r.get("qed", ""),
                "pains_ok": r.get("pains_ok", ""),
                "brenk_ok": r.get("brenk_ok", ""),
                "conformance_rate": r.get("conformance_rate", ""),
                "rotatable_between": r.get("rotatable_between", ""),
                "retro_feasibility": r.get("retro_feasibility", ""),
                "retro_n_cuts": r.get("retro_n_cuts", ""),
                "mw"      : r.get("mw", ""),
                "logp"    : r.get("logp", ""),
                "smiles"  : r.get("smiles", ""),
                "pair"    : r.get("pair", ""),
                "source"  : r.get("source", ""),
                "out_sdf" : r.get("out_sdf", ""),
            })
    print(f"\n  summary.csv 保存: {csv_path}")


# ──────────────────────────────────────────────────────────
# README 生成
# ──────────────────────────────────────────────────────────

def write_readme(selected: list[dict], out_dir: Path):
    lines = [
        "# Result_Best — 選抜分子一覧\n",
        "",
        "全ドッキング結果から **スコア上位 + Ligand Efficiency 上位** を選抜した最終候補。\n",
        "",
        "## 選抜基準",
        "| 基準 | 内容 |",
        "|------|------|",
        "| DrugLike | Lipinski Ro5 + Veber 充足 |",
        f"| スコアカットオフ | ≤ {SCORE_CUTOFF} kcal/mol |",
        f"| 統合上位 | Top {TOP_N} (スコア + LE の統合ランキング) |",
        "",
        "## コンフォメーション適合性指標",
        "| 指標 | 範囲 | 意味 |",
        "|------|------|------|",
        "| Conformance Rate | 0〜1 | 低エネルギーコンフォマーのうちファーマコフォア距離が許容範囲内の割合 (高い = 良好) |",
        "| Rotatable Between | 0〜 | ファーマコフォア間の回転可能結合数 (少ない = 剛直) |",
        "",
        "## 合成容易性指標",
        "| 指標 | 範囲 | 意味 |",
        "|------|------|------|",
        "| SA Score | 1〜10 | 合成容易性 (低い = 合成しやすい, ≤6 推奨) |",
        "| QED | 0〜1 | 薬剤らしさ (高い = 良好, ≥0.5 推奨) |",
        "| PAINS | OK/NG | アッセイ干渉構造なし = OK |",
        "| BRENK | OK/NG | 構造アラートなし = OK |",
        "",
        "## LE (Ligand Efficiency) の見方",
        "```",
        "LE = |smina Affinity| / 重原子数 (HAC)",
        "",
        "◎ LE ≥ 0.4  : 優秀 (フラグメント〜低分子の理想域)",
        "○ LE ≥ 0.3  : 良好 (低分子医薬品の目安)",
        "△ LE ≥ 0.2  : 許容範囲",
        "✕ LE < 0.2  : 非効率 (大きすぎ or 弱すぎ)",
        "```",
        "",
        "## 選抜分子一覧",
        "",
        f"| 順位 | スコア (kcal/mol) | HAC | LE | LE評価 | SA | QED | P | B | Conf | Rot | 分子名 |",
        f"|------|-----------------|-----|-----|--------|-----|------|---|---|------|-----|--------|",
    ]

    for rank, r in enumerate(selected, 1):
        import math
        le_str = f"{r['le']:.4f}" if r.get("le") is not None else "N/A"
        sa_v = r.get("sa_score", float("nan"))
        qed_v = r.get("qed", float("nan"))
        sa_str = f"{sa_v:.1f}" if sa_v is not None and not (isinstance(sa_v, float) and math.isnan(sa_v)) else "?"
        qed_str = f"{qed_v:.3f}" if qed_v is not None and not (isinstance(qed_v, float) and math.isnan(qed_v)) else "?"
        pains_str = "OK" if r.get("pains_ok", True) else "NG"
        brenk_str = "OK" if r.get("brenk_ok", True) else "NG"
        conf_v = r.get("conformance_rate")
        conf_str = f"{conf_v:.0%}" if conf_v is not None and not (isinstance(conf_v, float) and math.isnan(conf_v)) else "-"
        rot_v = r.get("rotatable_between")
        rot_str = str(rot_v) if rot_v is not None else "-"
        lines.append(
            f"| {rank} | {r['score']:.3f} | {r.get('hac','?')} | {le_str} "
            f"| {le_grade(r.get('le'))} | {sa_str} | {qed_str} | {pains_str} | {brenk_str} "
            f"| {conf_str} | {rot_str} | {r['name'][:40]} |"
        )

    lines += [
        "",
        "## ファイル構成",
        "```",
        "Result_Best/",
        "├── 01_score*_<name>.sdf        # ドッキングポーズ SDF (PyMOL/VMD で可視化可)",
        "├── ...                         # 02〜N 番まで同様",
        "├── retrosynthesis/             # 逆合成解析結果",
        "│   ├── 01_retrosynthesis_*.png # 合成スキーム図",
        "│   ├── ...                     # 各分子ごとの図",
        "│   └── retrosynthesis_detail.json  # 全詳細データ (JSON)",
        "├── summary.csv                 # 全選抜分子の詳細データ",
        "└── README.md                   # 本ファイル",
        "```",
        "",
        "## PyMOL での可視化",
        "```python",
        "load results/docking/receptor.pdb",
        "load Result_Best/01_score*.sdf",
        "# 受容体を表面表示",
        "show surface, receptor",
        "set transparency, 0.5, receptor",
        "```",
        "",
        "## 注意事項",
        "- smina (AutoDock Vina) スコアは**同一受容体での比較**に有効",
        "- ペプチド vs 低分子の直接比較は LE を参照すること",
        "- 実験的検証 (IC50, SPR, ITC 等) で確認が必要",
    ]

    readme_path = out_dir / "README.md"
    with open(readme_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"  README.md 保存: {readme_path}")


# ──────────────────────────────────────────────────────────
# 逆合成スキーム図生成
# ──────────────────────────────────────────────────────────

def _draw_mol_image(smiles: str, size: tuple = (300, 200),
                    highlight_color=None) -> "Image.Image":
    """SMILES から分子画像を生成"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        # 描画不能 → プレースホルダ
        img = Image.new("RGB", size, "white")
        d = ImageDraw.Draw(img)
        d.text((10, size[1] // 2 - 10), smiles[:30], fill="red")
        return img
    AllChem.Compute2DCoords(mol)
    drawer = Draw.MolDraw2DCairo(size[0], size[1])
    opts = drawer.drawOptions()
    opts.clearBackground = True
    if highlight_color:
        drawer.drawOptions().setHighlightColour(highlight_color)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    import io
    png_data = drawer.GetDrawingText()
    return Image.open(io.BytesIO(png_data))


def _try_get_font(size: int = 16):
    """利用可能なフォントを取得"""
    font_paths = [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Regular.ttf",
    ]
    for fp in font_paths:
        if os.path.exists(fp):
            try:
                return ImageFont.truetype(fp, size)
            except Exception:
                continue
    return ImageFont.load_default()


# 日本語テキスト → 英語への変換テーブル (フォント非対応時用)
_JA_TO_EN = {
    "C-C 単結合 (sp3)":                "C-C single bond (sp3)",
    "エーテル / エステル C-O":          "Ether/Ester C-O",
    "アミド C-N":                       "Amide C-N",
    "芳香族 C=C (ビアリール等)":        "Aromatic C=C (biaryl)",
    "C-N (アミン / スルホンアミド)":    "C-N (amine/sulfonamide)",
    "C-O (アルキル エーテル)":          "C-O (alkyl ether)",
    "C-N (ラクタム / 環状アミド)":      "C-N (lactam)",
    "N-N (ヒドラジン / ヒドラジド)":    "N-N (hydrazine/hydrazide)",
    "C-S (チオエーテル)":               "C-S (thioether)",
    "C=C (オレフィン メタセシス)":      "C=C (olefin metathesis)",
    "C-C (環外結合)":                   "C-C (exocyclic)",
    "芳香環 C-ヘテロ原子":             "Ar C-heteroatom",
    "C-N (複素環)":                     "C-N (heterocyclic)",
    "芳香環 C-C (Suzuki 等)":          "Ar C-C (Suzuki coupling)",
    "市販品レベル (SA≤3.0, そのまま購入可能な可能性大)": "Commercially available (SA<=3.0)",
    "市販品レベル (切断不要)":          "Commercially available (no decomposition needed)",
    "容易 (2 ステップ以内, 全フラグメント市販品レベル)":
        "Easy (<=2 steps, all fragments commercially available)",
    "容易 (1-2 ステップ, 市販フラグメントの組み合わせ)":
        "Easy (1-2 steps, commercially available fragments)",
    "比較的容易 (市販フラグメントの組み合わせ)":
        "Relatively easy (commercially available fragment combination)",
    "中程度 (3-4 ステップ, フラグメントは合成容易)":
        "Moderate (3-4 steps, fragments are easy to synthesize)",
    "中程度 (標準的な有機合成で到達可能)":
        "Moderate (achievable by standard organic synthesis)",
    "中程度 (やや複雑だが到達可能)":
        "Moderate (somewhat complex but achievable)",
    "やや困難 (特殊試薬 or 多段階合成)":
        "Somewhat difficult (specialized reagents or multi-step)",
    "やや困難 (特殊フラグメントを含む)":
        "Somewhat difficult (contains specialized fragments)",
    "困難 (5 ステップ以上の多段階合成)":
        "Difficult (5+ step multi-step synthesis)",
    "困難 (6 ステップ以上の多段階合成, 要専門家相談)":
        "Difficult (6+ steps, expert consultation recommended)",
    "合成可能 (AI 逆合成ルートあり)":
        "Synthesizable (AI retrosynthesis route found)",
    "解析不可":                         "Analysis not possible",
}


def _to_en(text: str) -> str:
    """日本語テキストが含まれていれば英語に変換"""
    if text in _JA_TO_EN:
        return _JA_TO_EN[text]
    # 部分一致も試す
    for ja, en in _JA_TO_EN.items():
        if ja in text:
            text = text.replace(ja, en)
    return text


# ──────────────────────────────────────────────────────────
# 前向き合成条件データベース
# ──────────────────────────────────────────────────────────

# BRICS 反応タイプ (日本語キー) → 前向き合成条件
FORWARD_SYNTHESIS = {
    "アミド C-N": {
        "name_en": "Amide Coupling",
        "reagents": "HATU, DIPEA (or EDC, HOBt, Et3N)",
        "solvent": "DMF or DCM",
        "temp": "0 °C → r.t., 2-12 h",
        "substrate": "R-COOH + R'-NH2",
        "workup": "Aqueous workup, column chromatography",
        "yield_range": "60-90%",
        "ref": "Standard peptide coupling",
    },
    "芳香環 C-C (Suzuki 等)": {
        "name_en": "Suzuki-Miyaura Cross-Coupling",
        "reagents": "Pd(PPh3)4 (5 mol%), K2CO3 (2 eq)",
        "solvent": "1,4-dioxane / H2O (3:1)",
        "temp": "80-100 °C, 12-24 h, N2 atm.",
        "substrate": "Ar-Br (or Ar-I) + Ar'-B(OH)2",
        "workup": "Filter through Celite, extract, column",
        "yield_range": "50-85%",
        "ref": "Miyaura, N.; Suzuki, A. Chem. Rev. 1995, 95, 2457",
    },
    "C-N (アミン / スルホンアミド)": {
        "name_en": "Reductive Amination",
        "reagents": "NaBH3CN (or NaBH(OAc)3), AcOH (cat.)",
        "solvent": "MeOH or DCE",
        "temp": "r.t., 2-6 h",
        "substrate": "R-CHO + R'-NH2",
        "workup": "NaHCO3 wash, extract, column",
        "yield_range": "50-80%",
        "ref": "Abdel-Magid, A. F. et al. J. Org. Chem. 1996, 61, 3849",
    },
    "C-N (複素環)": {
        "name_en": "N-Alkylation",
        "reagents": "K2CO3 or Cs2CO3 (2 eq)",
        "solvent": "DMF or MeCN",
        "temp": "60-80 °C, 6-12 h",
        "substrate": "Het-NH + R-X (X = Br, I, OTs)",
        "workup": "Dilute with H2O, extract with EtOAc, column",
        "yield_range": "50-75%",
        "ref": "Standard N-alkylation",
    },
    "エーテル / エステル C-O": {
        "name_en": "Esterification / Mitsunobu",
        "reagents": "DCC, DMAP (ester) or PPh3, DIAD (Mitsunobu)",
        "solvent": "DCM or THF",
        "temp": "0 °C → r.t., 4-12 h",
        "substrate": "R-COOH + R'-OH (ester) or R-OH + R'-OH (ether)",
        "workup": "Filter, wash, column",
        "yield_range": "60-85%",
        "ref": "Mitsunobu, O. Synthesis 1981, 1",
    },
    "C-O (アルキル エーテル)": {
        "name_en": "Williamson Ether Synthesis",
        "reagents": "NaH (1.2 eq), then R'-X (X = Br, I)",
        "solvent": "THF (dry)",
        "temp": "0 °C → r.t., 2-6 h",
        "substrate": "R-OH + R'-X",
        "workup": "Quench with H2O, extract, column",
        "yield_range": "60-80%",
        "ref": "Williamson, A. Phil. Mag. 1850, 37, 350",
    },
    "C-C 単結合 (sp3)": {
        "name_en": "Grignard Reaction / Alkylation",
        "reagents": "RMgBr (Grignard) or LDA + R'-X",
        "solvent": "THF (dry)",
        "temp": "-78 °C → r.t.",
        "substrate": "R-MgBr + R'-CHO (or ketone)",
        "workup": "NH4Cl quench, extract, column",
        "yield_range": "50-80%",
        "ref": "Standard Grignard / enolate alkylation",
    },
    "C-N (ラクタム / 環状アミド)": {
        "name_en": "Lactam Formation",
        "reagents": "T3P (propylphosphonic anhydride) or HATU",
        "solvent": "DMF, dilute conditions (0.01 M)",
        "temp": "r.t., 12-24 h",
        "substrate": "Amino acid (intramolecular cyclization)",
        "workup": "Dilute with H2O, extract, column",
        "yield_range": "30-60%",
        "ref": "Intramolecular amide cyclization",
    },
    "C=C (オレフィン メタセシス)": {
        "name_en": "Olefin Metathesis",
        "reagents": "Grubbs II catalyst (5-10 mol%)",
        "solvent": "DCM (degassed)",
        "temp": "40 °C, 6-24 h, N2 atm.",
        "substrate": "R-CH=CH2 + R'-CH=CH2",
        "workup": "Filter through SiO2, concentrate, column",
        "yield_range": "40-75%",
        "ref": "Grubbs, R. H. Angew. Chem. Int. Ed. 2006, 45, 3760",
    },
    "C-S (チオエーテル)": {
        "name_en": "Thioether Formation",
        "reagents": "K2CO3 (2 eq)",
        "solvent": "DMF",
        "temp": "r.t., 2-6 h",
        "substrate": "R-SH + R'-X (X = Br, I)",
        "workup": "Dilute with H2O, extract, column",
        "yield_range": "70-90%",
        "ref": "Standard S-alkylation",
    },
    "芳香環 C-ヘテロ原子": {
        "name_en": "Buchwald-Hartwig Amination",
        "reagents": "Pd2(dba)3 (2.5 mol%), XPhos, NaOtBu",
        "solvent": "toluene",
        "temp": "100 °C, 12-24 h, N2 atm.",
        "substrate": "Ar-Br + R-NH2",
        "workup": "Filter through Celite, column",
        "yield_range": "50-80%",
        "ref": "Surry, D. S.; Buchwald, S. L. Chem. Sci. 2011, 2, 27",
    },
    "C-C (環外結合)": {
        "name_en": "Heck Reaction / Negishi Coupling",
        "reagents": "Pd(OAc)2 (5 mol%), PPh3, Et3N",
        "solvent": "DMF",
        "temp": "80-120 °C, 12 h, N2 atm.",
        "substrate": "Ar-X + R-CH=CH2 (Heck) or R-ZnCl (Negishi)",
        "workup": "Filter, extract, column",
        "yield_range": "50-80%",
        "ref": "Heck, R. F. Org. React. 1982, 27, 345",
    },
    "N-N (ヒドラジン / ヒドラジド)": {
        "name_en": "Hydrazide Formation",
        "reagents": "R-COOH + NH2NH2·H2O, EDC",
        "solvent": "EtOH or DMF",
        "temp": "r.t. → 60 °C, 4-12 h",
        "substrate": "R-COOH + NH2-NHR'",
        "workup": "Concentrate, recrystallize or column",
        "yield_range": "60-85%",
        "ref": "Standard hydrazide synthesis",
    },
}

# BRICS type tuple → 日本語キーへのマッピング (retrosynthesis.py と対応)
_BRICS_TUPLE_TO_JA = {
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


def _get_synthesis_conditions(reaction_ja: str) -> dict | None:
    """BRICS 反応タイプ名(日本語)から前向き合成条件を取得"""
    if reaction_ja in FORWARD_SYNTHESIS:
        return FORWARD_SYNTHESIS[reaction_ja]
    # "BRICS type (X, Y)" 形式の場合
    import re
    m = re.match(r"BRICS type \((.+)\)", reaction_ja)
    if m:
        nums = tuple(sorted(int(x.strip()) for x in m.group(1).split(",")))
        # 各要素で個別検索
        for n in nums:
            ja_key = _BRICS_TUPLE_TO_JA.get((n,))
            if ja_key and ja_key in FORWARD_SYNTHESIS:
                return FORWARD_SYNTHESIS[ja_key]
    return None


def generate_synthesis_plan(retro_result: dict) -> list[dict]:
    """
    BRICS 逆合成結果から前向き合成プランを生成。

    Returns:
        [
          {
            "step": 1,
            "reaction_ja": "アミド C-N",
            "reaction_en": "Amide Coupling",
            "reagents": "HATU, DIPEA ...",
            "solvent": "DMF",
            "temp": "0°C → r.t., 2-12h",
            "substrate": "R-COOH + R'-NH2",
            "yield_range": "60-90%",
            "ref": "...",
            "fragments_used": ["CCC=O", "Nc1ccc2[nH]ccc2c1"],
          },
          ...
        ]
    """
    brics = retro_result.get("brics", {})
    reactions = brics.get("reactions", [])
    fragments = brics.get("fragments", [])
    n_cuts = brics.get("n_cuts", 0)

    if n_cuts == 0 or not reactions:
        return []

    # フラグメントをサイズ順 (大→小) にソート → 大きいのをコアとして使う
    sorted_frags = sorted(fragments, key=lambda f: f["n_atoms"], reverse=True)

    steps = []
    used_reactions = set()

    for i, rxn_ja in enumerate(reactions):
        if rxn_ja in used_reactions:
            continue
        used_reactions.add(rxn_ja)

        cond = _get_synthesis_conditions(rxn_ja)
        if cond is None:
            # 未知の反応タイプ → 汎用条件
            cond = {
                "name_en": _to_en(rxn_ja),
                "reagents": "See literature for specific conditions",
                "solvent": "Appropriate solvent",
                "temp": "Optimized conditions",
                "substrate": "Fragment A + Fragment B",
                "workup": "Standard workup",
                "yield_range": "Variable",
                "ref": "Reaction-type specific",
            }

        # どのフラグメントがこのステップで使われるか (近似的に割り当て)
        frags_for_step = []
        if i < len(sorted_frags):
            frags_for_step.append(sorted_frags[i]["smiles"])
        if i + 1 < len(sorted_frags):
            frags_for_step.append(sorted_frags[i + 1]["smiles"])

        steps.append({
            "step": len(steps) + 1,
            "reaction_ja": rxn_ja,
            "reaction_en": cond.get("name_en", _to_en(rxn_ja)),
            "reagents": cond["reagents"],
            "solvent": cond["solvent"],
            "temp": cond["temp"],
            "substrate": cond.get("substrate", ""),
            "workup": cond.get("workup", ""),
            "yield_range": cond.get("yield_range", ""),
            "ref": cond.get("ref", ""),
            "fragments_used": frags_for_step,
        })

    return steps


def draw_retrosynthesis_scheme(
    record: dict,
    retro_result: dict,
    out_path: Path,
    max_frags_per_row: int = 4,
) -> str | None:
    """
    1 分子の逆合成スキームを図示して PNG 保存。

    レイアウト:
      ┌─────────────────────────────┐
      │  Target Molecule (大きく)    │
      │  名前, SMILES, SA, QED      │
      ├─────────────────────────────┤
      │        ⇓ 逆合成              │
      │  反応タイプ: ...             │
      ├─────────────────────────────┤
      │  Fragment 1  Fragment 2  ...│
      │  SA=x.x      SA=x.x        │
      │  市販品       要合成         │
      ├─────────────────────────────┤
      │  総合判定: ...               │
      └─────────────────────────────┘
    """
    if not _HAS_DRAW:
        return None

    import math as _math

    smiles = record.get("smiles", "")
    name = record.get("name", "")
    brics = retro_result.get("brics", {})
    fragments = brics.get("fragments", [])
    reactions = brics.get("reactions", [])
    n_cuts = brics.get("n_cuts", 0)
    overall = retro_result.get("overall_feasibility", "?")

    font_large = _try_get_font(18)
    font_med   = _try_get_font(14)
    font_small = _try_get_font(12)

    # --- サイズ計算 ---
    mol_w, mol_h = 400, 280
    frag_w, frag_h = 250, 180
    n_frags = len(fragments)
    if n_frags == 0:
        n_frags = 1  # プレースホルダ
    n_cols = min(n_frags, max_frags_per_row)
    n_rows = _math.ceil(n_frags / max_frags_per_row)

    padding = 30
    frag_area_w = n_cols * (frag_w + 20) + padding * 2
    canvas_w = max(mol_w + padding * 2, frag_area_w, 600)
    # ヘッダ + 分子 + 矢印 + フラグメント + フッタ
    header_h = 50
    mol_section_h = mol_h + 80   # 分子 + テキスト
    arrow_h = 100                # 矢印 + 反応タイプ
    frag_section_h = n_rows * (frag_h + 70) + 20
    footer_h = 80
    canvas_h = header_h + mol_section_h + arrow_h + frag_section_h + footer_h

    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    draw = ImageDraw.Draw(canvas)

    y_cursor = padding

    # --- ヘッダ: タイトル ---
    title = f"Retrosynthesis Scheme: {name}"
    draw.text((padding, y_cursor), title, fill="black", font=font_large)
    y_cursor += header_h

    # --- ターゲット分子 ---
    target_img = _draw_mol_image(smiles, size=(mol_w, mol_h))
    mol_x = (canvas_w - mol_w) // 2
    canvas.paste(target_img, (mol_x, y_cursor))
    y_cursor += mol_h + 5

    # 分子情報テキスト
    sa_v = record.get("sa_score", "?")
    qed_v = record.get("qed", "?")
    sa_str = f"{sa_v:.1f}" if isinstance(sa_v, (int, float)) and not _math.isnan(sa_v) else "?"
    qed_str = f"{qed_v:.3f}" if isinstance(qed_v, (int, float)) and not _math.isnan(qed_v) else "?"
    info = f"SMILES: {smiles[:80]}{'...' if len(smiles) > 80 else ''}   |   SA Score: {sa_str}   |   QED: {qed_str}"
    draw.text((padding, y_cursor), info, fill="gray", font=font_small)
    y_cursor += 25
    score_str = f"Docking: {record.get('score', '?')} kcal/mol   |   LE: {record.get('le', '?')}"
    draw.text((padding, y_cursor), score_str, fill="gray", font=font_small)
    y_cursor += 30

    # --- 逆合成矢印 ---
    arrow_cx = canvas_w // 2
    arrow_top = y_cursor
    arrow_bottom = y_cursor + 40
    # 矢印本体
    draw.line([(arrow_cx, arrow_top), (arrow_cx, arrow_bottom)],
              fill="#2196F3", width=3)
    # 矢頭
    draw.polygon([
        (arrow_cx - 10, arrow_bottom),
        (arrow_cx + 10, arrow_bottom),
        (arrow_cx, arrow_bottom + 15),
    ], fill="#2196F3")
    y_cursor = arrow_bottom + 20

    # 反応タイプテキスト
    if reactions:
        en_reactions = [_to_en(r) for r in reactions]
        rxn_text = f"BRICS ({n_cuts} cuts): {', '.join(en_reactions[:3])}"
        if len(en_reactions) > 3:
            rxn_text += f"  (+{len(en_reactions) - 3} more)"
    else:
        rxn_text = f"BRICS: No decomposition (intact molecule)"
    draw.text((padding, y_cursor), rxn_text, fill="#1565C0", font=font_med)
    y_cursor += 35

    # --- フラグメント ---
    if not fragments:
        draw.text((padding, y_cursor), "No fragments (market-available compound)",
                  fill="#4CAF50", font=font_med)
        y_cursor += frag_h + 50
    else:
        for row_i in range(n_rows):
            start = row_i * max_frags_per_row
            end = min(start + max_frags_per_row, len(fragments))
            row_frags = fragments[start:end]
            row_w = len(row_frags) * (frag_w + 20) - 20
            x_start = (canvas_w - row_w) // 2

            for col_i, frag in enumerate(row_frags):
                fx = x_start + col_i * (frag_w + 20)
                fy = y_cursor

                # フラグメント枠
                buyable = frag.get("buyable", False)
                border_color = "#4CAF50" if buyable else "#FF9800"
                bg_color = "#F1F8E9" if buyable else "#FFF3E0"
                draw.rounded_rectangle(
                    [(fx - 5, fy - 5), (fx + frag_w + 5, fy + frag_h + 55)],
                    radius=8, fill=bg_color, outline=border_color, width=2)

                # フラグメント分子画像
                frag_img = _draw_mol_image(frag["smiles"], size=(frag_w, frag_h))
                canvas.paste(frag_img, (fx, fy))

                # フラグメント情報
                fsa = frag.get("sa_score", float("nan"))
                fsa_str = f"{fsa:.1f}" if not _math.isnan(fsa) else "?"
                buy_str = "Available" if buyable else "Needs synthesis"
                buy_color = "#2E7D32" if buyable else "#E65100"

                draw.text((fx + 5, fy + frag_h + 5),
                          f"SA={fsa_str}  ({frag['n_atoms']} atoms)",
                          fill="gray", font=font_small)
                draw.text((fx + 5, fy + frag_h + 22),
                          buy_str, fill=buy_color, font=font_small)
                # SMILES (短縮)
                frag_smi = frag["smiles"]
                if len(frag_smi) > 28:
                    frag_smi = frag_smi[:25] + "..."
                draw.text((fx + 5, fy + frag_h + 38),
                          frag_smi, fill="#666", font=font_small)

            y_cursor += frag_h + 75

    # --- フッタ: 総合判定 ---
    draw.line([(padding, y_cursor), (canvas_w - padding, y_cursor)],
              fill="#BDBDBD", width=1)
    y_cursor += 10

    # 色分け
    if "市販品" in overall or "容易" in overall:
        verdict_color = "#2E7D32"
    elif "中程度" in overall:
        verdict_color = "#F57F17"
    else:
        verdict_color = "#C62828"

    draw.text((padding, y_cursor),
              f"Feasibility: {_to_en(overall)}",
              fill=verdict_color, font=font_large)

    # 保存 (余白を切り詰め)
    actual_h = min(y_cursor + 40, canvas_h)
    canvas = canvas.crop((0, 0, canvas_w, actual_h))
    canvas.save(str(out_path), "PNG", quality=95)
    return str(out_path)


def _generate_synthesis_sequence(
    target_smiles: str,
    aizf_route: dict | None = None,
) -> list[dict]:
    """
    合成シーケンスを生成する。

    AiZynthFinder のルートがあればそれを使い、
    なければ従来の BRICS ベース分解にフォールバック。

    Returns:
        [{"step": 1,
          "substrate_smi": "...",
          "reagent_smi":   "..." or None,
          "product_smi":   "...",
          "reaction_en":   "...",  # 反応名 (AiZynthFinder の場合)
          "reagents":      "...",
          "solvent":       "...",
          "temp":          "...",
          "yield_range":   "...",
          "bond_types":    ("3","5") or None,  # BRICS の場合のみ
         }, ...]
    """
    # --- AiZynthFinder ルートが利用可能ならそちらを使用 ---
    if aizf_route and aizf_route.get("forward_steps"):
        return aizf_route["forward_steps"]

    # --- フォールバック: BRICS ベース分解 ---
    return _brics_based_synthesis_sequence(target_smiles)


def _brics_based_synthesis_sequence(target_smiles: str) -> list[dict]:
    """
    ターゲット分子のBRICS結合を段階的に形成し、
    各ステップの基質・試薬・中間体構造を生成する (従来の方式)。
    """
    from rdkit.Chem import BRICS, FragmentOnBonds

    mol = Chem.MolFromSmiles(target_smiles)
    if mol is None:
        return []

    brics_info = list(BRICS.FindBRICSBonds(mol))
    if not brics_info:
        return []

    bonds = []
    for (a1, a2), (t1, t2) in brics_info:
        if mol.GetBondBetweenAtoms(a1, a2):
            bonds.append({"a1": a1, "a2": a2, "types": (t1, t2)})
    if not bonds:
        return []

    # --- Helper: ターゲット分子の特定結合を切断し、フラグメントを返す ---
    def _break(bond_list):
        """bond_list の結合を切断 → [(atom_index_set, SMILES), ...] を返す"""
        if not bond_list:
            return [{"atoms": set(range(mol.GetNumAtoms())),
                     "smi": Chem.MolToSmiles(mol),
                     "size": mol.GetNumHeavyAtoms()}]

        bidxs = []
        for b in bond_list:
            bo = mol.GetBondBetweenAtoms(b["a1"], b["a2"])
            if bo:
                bidxs.append(bo.GetIdx())
        if not bidxs:
            return [{"atoms": set(range(mol.GetNumAtoms())),
                     "smi": Chem.MolToSmiles(mol),
                     "size": mol.GetNumHeavyAtoms()}]

        fragmented = FragmentOnBonds(mol, bidxs, addDummies=True)
        frag_mols = Chem.GetMolFrags(fragmented, asMols=True, sanitizeFrags=False)
        rwmol = Chem.RWMol(mol)
        for b in bond_list:
            if rwmol.GetBondBetweenAtoms(b["a1"], b["a2"]):
                rwmol.RemoveBond(b["a1"], b["a2"])
        frag_atom_tuples = Chem.GetMolFrags(rwmol)

        result = []
        for i, fmol in enumerate(frag_mols):
            edit = Chem.RWMol(fmol)
            for atom in edit.GetAtoms():
                if atom.GetAtomicNum() == 0:
                    atom.SetAtomicNum(1)
                    atom.SetIsotope(0)
                    atom.SetFormalCharge(0)
                    atom.SetNoImplicit(False)
            try:
                Chem.SanitizeMol(edit)
                edit = Chem.RemoveHs(edit)
                smi = Chem.MolToSmiles(edit)
            except Exception:
                smi = Chem.MolToSmiles(fmol)

            atoms_set = set(frag_atom_tuples[i]) if i < len(frag_atom_tuples) else set()
            result.append({"atoms": atoms_set, "smi": smi,
                           "size": len(atoms_set)})

        result.sort(key=lambda f: -f["size"])
        return result

    # --- 結合の合成順序を決定 (小さいフラグメントから先に結合) ---
    for b in bonds:
        frags = _break([b])
        b["min_frag"] = frags[-1]["size"] if frags else 0
    bonds.sort(key=lambda b: b["min_frag"])

    # --- 各ステップの中間体を生成 ---
    n = len(bonds)
    steps = []

    for k in range(n):
        frags_before = _break(bonds[k:])
        frags_after = _break(bonds[k + 1:])

        product = frags_after[0]

        a1, a2 = bonds[k]["a1"], bonds[k]["a2"]
        f_a1 = f_a2 = frags_before[0]
        for frag in frags_before:
            if a1 in frag["atoms"]:
                f_a1 = frag
            if a2 in frag["atoms"]:
                f_a2 = frag

        if f_a1["size"] >= f_a2["size"]:
            substrate, reagent = f_a1, f_a2
        else:
            substrate, reagent = f_a2, f_a1

        if substrate is reagent:
            reagent = None

        steps.append({
            "step":          k + 1,
            "substrate_smi": substrate["smi"],
            "reagent_smi":   reagent["smi"] if reagent else None,
            "product_smi":   product["smi"],
            "bond_types":    bonds[k]["types"],
        })

    return steps


def draw_forward_synthesis_scheme(
    record: dict,
    synth_seq: list[dict],
    out_path: Path,
) -> str | None:
    """
    Traditional organic chemistry synthesis scheme with intermediates:

    Step 1:
                       Amide Coupling
                     HATU, DIPEA / DMF, r.t.
      [SM_A] + [SM_B] ─────────────────────→ [Intermediate 1]
                          60-90%

    Step 2:
                       Suzuki Coupling
                     Pd(PPh3)4, K2CO3
      [Int 1] + [SM_C] ─────────────────────→ [Target]
                          50-85%
    """
    if not _HAS_DRAW:
        return None

    import math as _math, re as _re

    name = record.get("name", "")

    fonts = {
        "title":  _try_get_font(18),
        "step":   _try_get_font(13),
        "rxn":    _try_get_font(12),
        "cond":   _try_get_font(10),
        "plus":   _try_get_font(18),
        "label":  _try_get_font(9),
        "footer": _try_get_font(11),
    }

    # --- Layout constants ---
    MOL_W, MOL_H = 160, 120       # substrate / product image
    REAGENT_W, REAGENT_H = 120, 90  # reagent (smaller)
    ARROW_W = 200                  # arrow section width (shaft + text)
    ARROW_SHAFT = 160
    PLUS_W = 26
    GAP = 8
    PAD = 30
    TITLE_H = 50
    STEP_ROW_H = 200               # height per synthesis step row
    FOOTER_H = 35

    n_steps = len(synth_seq)

    # --- Pre-compute conditions for each step ---
    for s in synth_seq:
        bt = s.get("bond_types")
        # AiZynthFinder steps already have reaction_en set
        if s.get("reaction_en") and not bt:
            # Already filled by AiZynthFinder
            s.setdefault("reagents", "")
            s.setdefault("solvent", "")
            s.setdefault("temp", "")
            s.setdefault("yield_range", "")
            continue

        # BRICS-based: Map bond types → reaction conditions
        if bt:
            cond = None
            for t in bt:
                try:
                    t_int = int(t)
                except (ValueError, TypeError):
                    continue
                ja_key = _BRICS_TUPLE_TO_JA.get((t_int,))
                if ja_key:
                    cond = _get_synthesis_conditions(ja_key)
                    if cond:
                        s["reaction_en"] = cond.get("name_en", _to_en(ja_key))
                        s["reagents"]    = cond.get("reagents", "")
                        s["solvent"]     = cond.get("solvent", "")
                        s["temp"]        = cond.get("temp", "")
                        s["yield_range"] = cond.get("yield_range", "")
                        break
            if cond is None:
                s.setdefault("reaction_en", f"BRICS ({bt[0]},{bt[1]})")
                s.setdefault("reagents", "")
                s.setdefault("solvent", "")
                s.setdefault("temp", "")
                s.setdefault("yield_range", "")
        else:
            s.setdefault("reaction_en", "Transformation")
            s.setdefault("reagents", "")
            s.setdefault("solvent", "")
            s.setdefault("temp", "")
            s.setdefault("yield_range", "")

    # --- Canvas size ---
    # Each step row: [substrate] (+ [reagent]) [arrow] [product]
    row_w = MOL_W + PLUS_W + REAGENT_W + ARROW_W + MOL_W + GAP * 5
    canvas_w = max(row_w + PAD * 2, 800)
    canvas_h = TITLE_H + n_steps * STEP_ROW_H + FOOTER_H + PAD * 2

    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")
    draw = ImageDraw.Draw(canvas)

    # --- Title ---
    draw.text((PAD, PAD), f"Synthesis: {name}", fill="black", font=fonts["title"])

    y = PAD + TITLE_H

    # --- Draw helper: centered text above arrow ---
    def _ctxt(cx, ty, text, font, color):
        bb = draw.textbbox((0, 0), text, font=font)
        tw = bb[2] - bb[0]
        draw.text((cx - tw // 2, ty), text, fill=color, font=font)

    # --- Draw each step ---
    for s in synth_seq:
        step_n = s["step"]
        has_reagent = s["reagent_smi"] is not None
        is_last = (step_n == n_steps)

        # Step label
        draw.text((PAD, y), f"Step {step_n}", fill="#888", font=fonts["step"])
        y += 20

        # Vertical center for arrow within this step
        arrow_cy = y + (STEP_ROW_H - 40) // 2
        x = PAD + 10

        # --- Substrate ---
        sub_top = arrow_cy - MOL_H // 2
        sub_img = _draw_mol_image(s["substrate_smi"], size=(MOL_W, MOL_H))
        canvas.paste(sub_img, (x, sub_top))
        x += MOL_W + GAP

        # --- "+" and reagent ---
        if has_reagent:
            draw.text((x + 2, arrow_cy - 10), "+",
                      fill="#333", font=fonts["plus"])
            x += PLUS_W + GAP

            rea_top = arrow_cy - REAGENT_H // 2
            rea_img = _draw_mol_image(s["reagent_smi"],
                                      size=(REAGENT_W, REAGENT_H))
            canvas.paste(rea_img, (x, rea_top))
            x += REAGENT_W + GAP
        else:
            x += GAP

        # --- Arrow with conditions ---
        acx = x + ARROW_W // 2   # center of arrow text area

        # Reaction name
        rxn = s.get("reaction_en", "")
        if len(rxn) > 36:
            rxn = rxn[:33] + "..."
        _ctxt(acx, arrow_cy - 42, rxn, fonts["rxn"], "#C62828")

        # Reagents
        reagents = s.get("reagents", "")
        if len(reagents) > 38:
            reagents = reagents[:35] + "..."
        _ctxt(acx, arrow_cy - 26, reagents, fonts["cond"], "#333")

        # Solvent, temp
        sol_temp = f"{s.get('solvent', '')}, {s.get('temp', '')}"
        if len(sol_temp) > 38:
            sol_temp = sol_temp[:35] + "..."
        _ctxt(acx, arrow_cy - 12, sol_temp, fonts["cond"], "#555")

        # Arrow shaft
        ax1 = x + (ARROW_W - ARROW_SHAFT) // 2
        ax2 = ax1 + ARROW_SHAFT
        draw.line([(ax1, arrow_cy), (ax2, arrow_cy)], fill="black", width=2)

        # Arrowhead
        draw.polygon([
            (ax2, arrow_cy),
            (ax2 - 12, arrow_cy - 6),
            (ax2 - 12, arrow_cy + 6),
        ], fill="black")

        # Yield
        yld = s.get("yield_range", "")
        if yld:
            _ctxt(acx, arrow_cy + 8, yld, fonts["label"], "#888")

        x += ARROW_W + GAP

        # --- Product (intermediate or target) ---
        prod_top = arrow_cy - MOL_H // 2
        prod_img = _draw_mol_image(s["product_smi"], size=(MOL_W, MOL_H))

        # Blue border for final target
        if is_last:
            draw.rounded_rectangle(
                [(x - 4, prod_top - 4),
                 (x + MOL_W + 4, prod_top + MOL_H + 4)],
                radius=6, outline="#1565C0", width=2)

        canvas.paste(prod_img, (x, prod_top))

        # Label under product
        if is_last:
            sa_v = record.get("sa_score", "?")
            sa_s = (f"SA={sa_v:.1f}"
                    if isinstance(sa_v, (int, float))
                    and not _math.isnan(sa_v) else "")
            score_v = record.get("score", "?")
            draw.text((x, prod_top + MOL_H + 8),
                      f"Target  Score: {score_v}  {sa_s}",
                      fill="#1565C0", font=fonts["label"])

        y += STEP_ROW_H - 20

    # --- Footer: overall yield estimate ---
    total_low, total_high = 1.0, 1.0
    for s in synth_seq:
        nums = _re.findall(r"(\d+)", s.get("yield_range", ""))
        if len(nums) >= 2:
            total_low *= int(nums[0]) / 100
            total_high *= int(nums[1]) / 100
    est = (f"Overall: {n_steps} steps, "
           f"est. {total_low * 100:.0f}-{total_high * 100:.0f}% yield")
    draw.text((PAD, y), est, fill="#666", font=fonts["footer"])
    y += 25

    # --- Crop and save ---
    actual_h = min(y + PAD, canvas_h)
    canvas = canvas.crop((0, 0, canvas_w, actual_h))
    canvas.save(str(out_path), "PNG", quality=95)
    return str(out_path)


def run_detailed_retrosynthesis(
    selected: list[dict],
    out_dir: Path,
    use_aizynthfinder: bool = True,
    aizynthfinder_config: str = "",
):
    """選抜分子に対して詳細な逆合成解析を実行し、スキーム図を生成"""
    if not _HAS_RETRO:
        print("  [警告] retrosynthesis モジュールが利用不可")
        return

    retro_dir = out_dir / "retrosynthesis"
    retro_dir.mkdir(exist_ok=True)

    # AiZynthFinder の利用可否
    do_aizf = use_aizynthfinder and _HAS_AIZF
    if use_aizynthfinder and not _HAS_AIZF:
        print("  [警告] AiZynthFinder が利用不可 → BRICS フォールバック")

    print("\n  ④ 逆合成解析 + スキーム図生成...")
    if do_aizf:
        print("     (AiZynthFinder AI ベース逆合成を使用)")

    all_retro = []
    for rank, r in enumerate(selected, 1):
        smiles = r.get("smiles", "")
        name = r.get("name", "")
        if not smiles:
            continue

        print(f"     [{rank}] {name}...")
        retro = _retro_analyze(
            smiles, name=name, verbose=False,
            use_aizynthfinder=do_aizf,
            aizynthfinder_config=aizynthfinder_config,
        )
        all_retro.append(retro)

        # AiZynthFinder のベストルートを取得
        aizf_route = None
        if do_aizf and retro.get("aizynthfinder"):
            aizf_data = retro["aizynthfinder"]
            if aizf_data.get("routes"):
                aizf_route = aizf_data["routes"][0]
                n_steps = aizf_route.get("n_steps", "?")
                score = aizf_route.get("score", "?")
                print(f"         AiZynthFinder: {n_steps} steps, score={score}")

        # 逆合成スキーム図
        if _HAS_DRAW:
            img_name = f"{rank:02d}_retrosynthesis_{safe_filename(name)}.png"
            img_path = retro_dir / img_name
            result_path = draw_retrosynthesis_scheme(r, retro, img_path)
            if result_path:
                print(f"         → {img_name}")
                r["retro_scheme_png"] = result_path

        # 前向き合成スキーム図 (中間体付き)
        synth_steps = generate_synthesis_plan(retro)   # for HTML text display
        r["synth_steps"] = synth_steps
        if _HAS_DRAW:
            synth_seq = _generate_synthesis_sequence(
                smiles, aizf_route=aizf_route)
            if synth_seq:
                synth_name = f"{rank:02d}_synthesis_{safe_filename(name)}.png"
                synth_path = retro_dir / synth_name
                synth_result = draw_forward_synthesis_scheme(
                    r, synth_seq, synth_path)
                if synth_result:
                    src = "AiZynthFinder" if aizf_route else "BRICS"
                    print(f"         → {synth_name} ({src})")
                    r["synth_scheme_png"] = synth_result

        # 詳細情報を record に格納
        r["retro_detail"] = retro

    # 逆合成サマリー JSON 保存
    import json
    retro_json = []
    for retro in all_retro:
        entry = {
            "name":   retro["name"],
            "smiles": retro["smiles"],
            "overall_feasibility": retro["overall_feasibility"],
            "brics": {
                "n_cuts":      retro["brics"]["n_cuts"],
                "feasibility": retro["brics"]["feasibility"],
                "reactions":   retro["brics"]["reactions"],
                "fragments": [
                    {"smiles": f["smiles"], "sa_score": f["sa_score"],
                     "n_atoms": f["n_atoms"], "buyable": f["buyable"]}
                    for f in retro["brics"]["fragments"]
                ],
            },
            "recap": {
                "tree_depth": retro["recap"]["tree_depth"],
                "n_leaves":   len(retro["recap"]["all_leaves"]),
            },
        }
        retro_json.append(entry)

    json_path = retro_dir / "retrosynthesis_detail.json"
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(retro_json, f, indent=2, ensure_ascii=False)
    print(f"     → retrosynthesis_detail.json 保存")

    # HTML レポート生成
    html_path = _write_retrosynthesis_html(selected, all_retro, retro_dir)
    if html_path:
        print(f"     → {Path(html_path).name} 保存")

    return all_retro


def _img_to_base64(img_path: str) -> str:
    """画像を base64 エンコードしてインライン埋め込み用文字列を返す"""
    import base64
    with open(img_path, "rb") as f:
        data = base64.b64encode(f.read()).decode("ascii")
    return f"data:image/png;base64,{data}"


def _mol_svg_inline(smiles: str, w: int = 250, h: int = 180) -> str:
    """SMILES → SVG 文字列 (HTML インライン埋め込み用)"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return f'<span style="color:red">{smiles}</span>'
    AllChem.Compute2DCoords(mol)
    drawer = Draw.MolDraw2DSVG(w, h)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def _feasibility_badge(text: str) -> str:
    """合成可能性テキスト → HTML バッジ"""
    en = _to_en(text)
    if "市販品" in text or "Commercially" in en:
        color, bg = "#1B5E20", "#E8F5E9"
    elif "容易" in text or "Easy" in en:
        color, bg = "#2E7D32", "#F1F8E9"
    elif "比較的" in text or "Relatively" in en:
        color, bg = "#558B2F", "#F1F8E9"
    elif "中程度" in text or "Moderate" in en:
        color, bg = "#F57F17", "#FFF8E1"
    else:
        color, bg = "#C62828", "#FFEBEE"
    return (f'<span style="display:inline-block;padding:4px 12px;border-radius:4px;'
            f'background:{bg};color:{color};font-weight:600;'
            f'border:1px solid {color}">{en}</span>')


def _write_retrosynthesis_html(
    selected: list[dict],
    all_retro: list[dict],
    retro_dir: Path,
) -> str | None:
    """逆合成結果を見やすい HTML レポートとして出力"""
    import math as _math

    html_parts = []

    # --- CSS + ヘッダ ---
    html_parts.append("""<!DOCTYPE html>
<html lang="ja">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Retrosynthesis Report</title>
<style>
  * { box-sizing: border-box; margin: 0; padding: 0; }
  body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto,
         "Helvetica Neue", Arial, sans-serif;
         background: #f5f5f5; color: #333; line-height: 1.6; }
  .container { max-width: 1100px; margin: 0 auto; padding: 20px; }
  h1 { text-align: center; margin: 30px 0 10px; color: #1565C0; }
  .subtitle { text-align: center; color: #666; margin-bottom: 30px; }
  .mol-card { background: #fff; border-radius: 12px; padding: 24px;
              margin-bottom: 30px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); }
  .mol-header { display: flex; align-items: center; gap: 12px;
                margin-bottom: 16px; flex-wrap: wrap; }
  .mol-rank { display: inline-flex; align-items: center; justify-content: center;
              width: 36px; height: 36px; border-radius: 50%;
              background: #1565C0; color: #fff; font-weight: bold; font-size: 16px; }
  .mol-name { font-size: 20px; font-weight: 600; }
  .mol-info { display: grid; grid-template-columns: repeat(auto-fit, minmax(140px, 1fr));
              gap: 8px 16px; margin: 12px 0; padding: 12px;
              background: #FAFAFA; border-radius: 8px; font-size: 14px; }
  .mol-info dt { color: #888; font-size: 12px; }
  .mol-info dd { font-weight: 600; margin: 0; }
  .smiles-box { font-family: "Courier New", monospace; font-size: 12px;
                background: #f0f4f8; padding: 8px 12px; border-radius: 6px;
                word-break: break-all; color: #555; margin: 8px 0; }
  .scheme-img { text-align: center; margin: 20px 0; }
  .scheme-img img { max-width: 100%; border: 1px solid #e0e0e0; border-radius: 8px; }
  .section-title { font-size: 16px; font-weight: 600; color: #1565C0;
                   margin: 20px 0 10px; padding-bottom: 4px;
                   border-bottom: 2px solid #E3F2FD; }
  .frag-grid { display: flex; flex-wrap: wrap; gap: 16px; margin: 12px 0; }
  .frag-card { border-radius: 8px; padding: 12px; min-width: 200px; flex: 1;
               max-width: 280px; text-align: center; }
  .frag-card.buyable { background: #F1F8E9; border: 2px solid #4CAF50; }
  .frag-card.needs-synth { background: #FFF3E0; border: 2px solid #FF9800; }
  .frag-card svg { max-width: 100%; height: auto; }
  .frag-label { font-size: 12px; margin-top: 8px; }
  .frag-smiles { font-family: monospace; font-size: 11px; color: #777;
                 word-break: break-all; }
  .tag { display: inline-block; padding: 2px 8px; border-radius: 4px;
         font-size: 12px; font-weight: 600; }
  .tag-green { background: #E8F5E9; color: #2E7D32; }
  .tag-orange { background: #FFF3E0; color: #E65100; }
  .reactions-list { list-style: none; padding: 0; }
  .reactions-list li { padding: 4px 0; font-size: 14px; }
  .reactions-list li::before { content: "\\2192 "; color: #1565C0; font-weight: bold; }
  .feasibility-box { margin: 16px 0; padding: 12px; border-radius: 8px;
                     background: #FAFAFA; }
  .summary-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
  .summary-table th, .summary-table td { padding: 10px 12px; text-align: left;
                                          border-bottom: 1px solid #e0e0e0; }
  .summary-table th { background: #E3F2FD; color: #1565C0; font-size: 13px; }
  .summary-table tr:hover { background: #f5f5f5; }
  .arrow-down { text-align: center; font-size: 36px; color: #1565C0; margin: 8px 0; }
</style>
</head>
<body>
<div class="container">
<h1>Retrosynthesis Report</h1>
<p class="subtitle">Top 5 candidates &mdash; BRICS decomposition &amp; synthesis feasibility</p>
""")

    # --- サマリーテーブル ---
    html_parts.append("""
<table class="summary-table">
<tr>
  <th>#</th><th>Name</th><th>Score (kcal/mol)</th><th>LE</th>
  <th>SA</th><th>QED</th><th>BRICS cuts</th><th>Feasibility</th>
</tr>
""")
    for rank, (rec, retro) in enumerate(zip(selected, all_retro), 1):
        sa = rec.get("sa_score", "")
        sa_s = f"{sa:.1f}" if isinstance(sa, (int, float)) and not _math.isnan(sa) else "?"
        qed = rec.get("qed", "")
        qed_s = f"{qed:.3f}" if isinstance(qed, (int, float)) and not _math.isnan(qed) else "?"
        le = rec.get("le", "")
        le_s = f"{le:.4f}" if le is not None else "N/A"
        badge = _feasibility_badge(retro["overall_feasibility"])
        html_parts.append(
            f'<tr><td><strong>{rank}</strong></td>'
            f'<td><a href="#mol-{rank}" style="color:#1565C0">{retro["name"]}</a></td>'
            f'<td>{rec.get("score","")}</td>'
            f'<td>{le_s}</td>'
            f'<td>{sa_s}</td><td>{qed_s}</td>'
            f'<td>{retro["brics"]["n_cuts"]}</td>'
            f'<td>{badge}</td></tr>\n'
        )
    html_parts.append("</table>\n<hr style='margin:30px 0'>\n")

    # --- 各分子の詳細カード ---
    for rank, (rec, retro) in enumerate(zip(selected, all_retro), 1):
        sa = rec.get("sa_score", "")
        sa_s = f"{sa:.1f}" if isinstance(sa, (int, float)) and not _math.isnan(sa) else "?"
        qed = rec.get("qed", "")
        qed_s = f"{qed:.3f}" if isinstance(qed, (int, float)) and not _math.isnan(qed) else "?"
        le = rec.get("le", "")
        le_s = f"{le:.4f}" if le is not None else "N/A"
        le_g = le_grade(le)

        brics = retro["brics"]
        recap = retro["recap"]
        frags = brics.get("fragments", [])
        reactions = brics.get("reactions", [])

        html_parts.append(f"""
<div class="mol-card" id="mol-{rank}">
  <div class="mol-header">
    <span class="mol-rank">{rank}</span>
    <span class="mol-name">{retro['name']}</span>
  </div>

  <div class="smiles-box">{retro['smiles']}</div>

  <div class="mol-info">
    <div><dt>Docking Score</dt><dd>{rec.get('score','')} kcal/mol</dd></div>
    <div><dt>Ligand Efficiency</dt><dd>{le_s} ({le_g})</dd></div>
    <div><dt>SA Score</dt><dd>{sa_s}</dd></div>
    <div><dt>QED</dt><dd>{qed_s}</dd></div>
    <div><dt>Heavy Atoms</dt><dd>{rec.get('hac','?')}</dd></div>
    <div><dt>PAINS</dt><dd>{'OK' if rec.get('pains_ok', True) else 'NG'}</dd></div>
    <div><dt>BRENK</dt><dd>{'OK' if rec.get('brenk_ok', True) else 'NG'}</dd></div>
    <div><dt>Source</dt><dd>{rec.get('source','')}</dd></div>
  </div>
""")

        # 前向き合成スキーム図 (メインで表示)
        synth_path = rec.get("synth_scheme_png")
        if synth_path and os.path.exists(synth_path):
            b64_synth = _img_to_base64(synth_path)
            html_parts.append(f"""
  <div class="section-title" style="color:#E65100; border-bottom-color:#FFF3E0; font-size:18px">
    Synthesis Scheme
  </div>
  <div class="scheme-img">
    <img src="{b64_synth}" alt="Synthesis scheme for {retro['name']}">
  </div>
""")

        # 逆合成スキーム図 (折りたたみ)
        scheme_path = rec.get("retro_scheme_png")
        if scheme_path and os.path.exists(scheme_path):
            b64_retro = _img_to_base64(scheme_path)
            html_parts.append(f"""
  <details style="margin:10px 0">
    <summary style="cursor:pointer; color:#1565C0; font-weight:600; font-size:14px">
      Show Retrosynthetic Decomposition
    </summary>
    <div class="scheme-img" style="margin-top:10px">
      <img src="{b64_retro}" alt="Retrosynthesis scheme for {retro['name']}">
    </div>
  </details>
""")

        # BRICS 反応タイプ
        html_parts.append(f"""
  <div class="section-title">BRICS Decomposition ({brics['n_cuts']} cuts)</div>
""")
        if reactions:
            html_parts.append('  <ul class="reactions-list">\n')
            for rxn in reactions:
                html_parts.append(f'    <li>{_to_en(rxn)}</li>\n')
            html_parts.append('  </ul>\n')
        else:
            html_parts.append('  <p style="color:#666">No BRICS decomposition (intact molecule)</p>\n')

        # フラグメント
        if frags:
            html_parts.append("""
  <div class="section-title">Fragments</div>
  <div class="frag-grid">
""")
            for fi, frag in enumerate(frags, 1):
                buyable = frag.get("buyable", False)
                css_class = "buyable" if buyable else "needs-synth"
                fsa = frag.get("sa_score", float("nan"))
                fsa_s = f"{fsa:.1f}" if not _math.isnan(fsa) else "?"
                tag_class = "tag-green" if buyable else "tag-orange"
                tag_text = "Commercially available" if buyable else "Needs synthesis"
                svg = _mol_svg_inline(frag["smiles"], 200, 140)

                html_parts.append(f"""
    <div class="frag-card {css_class}">
      {svg}
      <div class="frag-label">
        <strong>Fragment {fi}</strong> &middot; {frag['n_atoms']} atoms &middot; SA = {fsa_s}
      </div>
      <div><span class="tag {tag_class}">{tag_text}</span></div>
      <div class="frag-smiles">{frag['smiles']}</div>
    </div>
""")
            html_parts.append("  </div>\n")

        # ── 前向き合成スキーム ──
        synth_steps = generate_synthesis_plan(retro)
        if synth_steps:
            html_parts.append("""
  <div class="section-title" style="color:#E65100; border-bottom-color:#FFF3E0">
    Proposed Synthesis Route
  </div>
""")
            for step in synth_steps:
                html_parts.append(f"""
  <div style="background:#FFFDE7; border-left:4px solid #F9A825;
              border-radius:0 8px 8px 0; padding:14px 18px; margin:10px 0">
    <div style="display:flex; align-items:center; gap:10px; margin-bottom:8px">
      <span style="display:inline-flex; align-items:center; justify-content:center;
                   min-width:28px; height:28px; border-radius:50%;
                   background:#F9A825; color:#fff; font-weight:bold; font-size:14px">
        {step['step']}
      </span>
      <strong style="font-size:16px">{step['reaction_en']}</strong>
    </div>
    <table style="font-size:13px; border-collapse:collapse; width:100%">
      <tr><td style="padding:3px 8px; color:#888; width:100px">Substrate</td>
          <td style="padding:3px 8px; font-family:monospace">{step['substrate']}</td></tr>
      <tr><td style="padding:3px 8px; color:#888">Reagents</td>
          <td style="padding:3px 8px"><strong>{step['reagents']}</strong></td></tr>
      <tr><td style="padding:3px 8px; color:#888">Solvent</td>
          <td style="padding:3px 8px">{step['solvent']}</td></tr>
      <tr><td style="padding:3px 8px; color:#888">Conditions</td>
          <td style="padding:3px 8px">{step['temp']}</td></tr>
      <tr><td style="padding:3px 8px; color:#888">Workup</td>
          <td style="padding:3px 8px">{step.get('workup','')}</td></tr>
      <tr><td style="padding:3px 8px; color:#888">Expected yield</td>
          <td style="padding:3px 8px">{step['yield_range']}</td></tr>
      <tr><td style="padding:3px 8px; color:#888">Reference</td>
          <td style="padding:3px 8px; font-style:italic; font-size:12px">{step['ref']}</td></tr>
    </table>
  </div>
""")
            # 合成全体の見積もり
            n_steps = len(synth_steps)
            # 各ステップの yield 中央値で推定
            import re as _re
            total_yield_low, total_yield_high = 1.0, 1.0
            for step in synth_steps:
                yr = step.get("yield_range", "")
                nums = _re.findall(r"(\d+)", yr)
                if len(nums) >= 2:
                    total_yield_low *= int(nums[0]) / 100
                    total_yield_high *= int(nums[1]) / 100
                elif len(nums) == 1:
                    total_yield_low *= int(nums[0]) / 100
                    total_yield_high *= int(nums[0]) / 100

            html_parts.append(f"""
  <div style="background:#E8EAF6; border-radius:8px; padding:12px 18px;
              margin:12px 0; font-size:14px">
    <strong>Overall estimate:</strong>
    {n_steps} step{'s' if n_steps > 1 else ''} &middot;
    estimated overall yield: {total_yield_low*100:.0f}-{total_yield_high*100:.0f}%
  </div>
""")
        else:
            html_parts.append("""
  <div class="section-title" style="color:#E65100; border-bottom-color:#FFF3E0">
    Proposed Synthesis Route
  </div>
  <p style="font-size:14px; color:#555; padding:8px 0">
    No BRICS decomposition &mdash; this compound may need to be obtained
    commercially or synthesized via a custom total synthesis route.
  </p>
""")

        # 総合判定
        badge = _feasibility_badge(retro["overall_feasibility"])
        html_parts.append(f"""
  <div class="feasibility-box">
    <strong>Overall Feasibility:</strong> {badge}
  </div>
</div>
""")

    # --- フッタ ---
    html_parts.append("""
<footer style="text-align:center; color:#999; font-size:12px; margin:40px 0 20px">
  Generated by Peptide-to-Small-Molecule Pipeline &mdash; Retrosynthesis Module (BRICS / RECAP)
</footer>
</div>
</body>
</html>
""")

    html_path = retro_dir / "retrosynthesis_report.html"
    with open(html_path, "w", encoding="utf-8") as f:
        f.write("".join(html_parts))

    return str(html_path)


# ──────────────────────────────────────────────────────────
# メイン
# ──────────────────────────────────────────────────────────

def main():
    print("\n" + "=" * 60)
    print("  Result_Best 選抜スクリプト")
    print("=" * 60)

    # 出力ディレクトリ初期化 (古いファイルを削除)
    if OUT_DIR.exists():
        for old_file in OUT_DIR.glob("*.sdf"):
            old_file.unlink()
        for old_file in OUT_DIR.glob("*.csv"):
            old_file.unlink()
        for old_file in OUT_DIR.glob("*.md"):
            old_file.unlink()
    OUT_DIR.mkdir(exist_ok=True)
    print(f"  出力先: {OUT_DIR.resolve()}\n")

    # 全結果読み込み
    print("  ① 全ドッキング結果を読み込み中...")
    records = load_all_results()
    print(f"     合計 {len(records)} 件 (ドッキングスコアあり)")

    # 選抜
    print(f"\n  ② 選抜 (DrugLike + スコア ≤ {SCORE_CUTOFF} + 統合 Top{TOP_N})")
    selected = select_best(records)
    print(f"     → {len(selected)} 件 選抜\n")

    if not selected:
        print("  選抜分子なし。基準を緩めるか results/ を確認してください。")
        return

    # SDF コピー
    print("  ③ SDF ファイルをコピー中...")
    selected = copy_sdfs(selected, OUT_DIR)

    # 逆合成解析 + スキーム図 (AiZynthFinder があれば自動使用)
    run_detailed_retrosynthesis(selected, OUT_DIR, use_aizynthfinder=_HAS_AIZF)

    # サマリー CSV
    write_summary_csv(selected, OUT_DIR)

    # README
    write_readme(selected, OUT_DIR)

    # ─── 最終サマリー表示 ───
    import math
    print()
    print("=" * 160)
    print(f"  {'順位':<5} {'スコア':>9}  {'HAC':>5}  {'LE':>7}  {'LE評価':<10}  {'SA':>5}  {'QED':>6}  {'P':>2} {'B':>2}  {'Conf':>5} {'Rot':>3}  {'逆合成':>12}  分子名")
    print("=" * 160)
    for rank, r in enumerate(selected, 1):
        le_str  = f"{r['le']:.4f}" if r.get("le") is not None else "  N/A "
        grade   = le_grade(r.get("le"))
        if "STEP3" in r["source"]:
            src_tag = "STEP3"
        elif "STEP7b" in r["source"]:
            src_tag = "STEP7b"
        else:
            src_tag = "STEP7"
        sa_v = r.get("sa_score")
        qed_v = r.get("qed")
        sa_str = f"{sa_v:.1f}" if sa_v is not None and not (isinstance(sa_v, float) and math.isnan(sa_v)) else "  ? "
        qed_str = f"{qed_v:.3f}" if qed_v is not None and not (isinstance(qed_v, float) and math.isnan(qed_v)) else "  ?  "
        p_str = "OK" if r.get("pains_ok", True) else "NG"
        b_str = "OK" if r.get("brenk_ok", True) else "NG"
        conf_v = r.get("conformance_rate")
        conf_str = f"{conf_v:.0%}" if conf_v is not None and not (isinstance(conf_v, float) and math.isnan(conf_v)) else "  -  "
        rot_v = r.get("rotatable_between")
        rot_str = f"{rot_v:>3}" if rot_v is not None else "  -"
        retro = r.get("retro_feasibility", "?")
        # 短縮表示
        if "市販品" in retro:
            retro_short = "市販品"
        elif "容易" in retro and "比較的" not in retro:
            retro_short = "容易"
        elif "比較的" in retro:
            retro_short = "比較的容易"
        elif "中程度" in retro:
            retro_short = "中程度"
        elif "困難" in retro and "やや" in retro:
            retro_short = "やや困難"
        elif "困難" in retro:
            retro_short = "困難"
        else:
            retro_short = retro[:8]
        n_cuts = r.get("retro_n_cuts", "?")
        print(f"  {rank:<5} {r['score']:>9.3f}  {str(r.get('hac','?')):>5}  {le_str:>7}  "
              f"{grade:<10}  {sa_str:>5}  {qed_str:>6}  {p_str:>2} {b_str:>2}  "
              f"{conf_str:>5} {rot_str:>3}  "
              f"{retro_short:<10}({n_cuts})  [{src_tag}] {r['name'][:35]}")
    print("=" * 160)
    print(f"\n  完了! Result_Best/ に {len(selected)} 件のSDF + summary.csv + README.md を保存しました")


if __name__ == "__main__":
    main()
