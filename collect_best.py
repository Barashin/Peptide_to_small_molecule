#!/usr/bin/env python3
"""
collect_best.py
===============
全ドッキング結果を集約し、結合力・Ligand Efficiency の上位分子を
Result_Best/ ディレクトリにまとめる。

選抜基準:
  1. dock_score あり + ドッキング SDF が実在
  2. DrugLike = True (bridge 結果のみ判定可; peptidomimetics は常に含める)
  3. smina スコア上位 + LE 上位の 2 軸でランキング
  4. Top_Score 10 件 + Top_LE 5 件 (重複除去) を Result_Best/ にコピー

出力:
  Result_Best/
  ├── <rank>_<score>_<name>.sdf   # ドッキングポーズ SDF
  ├── summary.csv                 # 全選抜分子の詳細
  └── README.md                   # 見方と考察
"""

import os
import csv
import shutil
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import Descriptors


# ──────────────────────────────────────────────────────────
# 設定
# ──────────────────────────────────────────────────────────

RESULTS_DIR  = Path("results")
OUT_DIR      = Path("Result_Best")
TOP_SCORE_N  = 10   # スコア上位 N 件
TOP_LE_N     = 5    # LE 上位 N 件 (スコア上位と重複除去)
SCORE_CUTOFF = -5.0 # これより弱い結合は除外


# ──────────────────────────────────────────────────────────
# ユーティリティ
# ──────────────────────────────────────────────────────────

def calc_hac(smiles: str) -> int | None:
    if not smiles:
        return None
    mol = Chem.MolFromSmiles(smiles)
    return mol.GetNumAtoms() if mol else None


def calc_le(score: float | None, hac: int | None) -> float | None:
    if score is None or hac is None or hac == 0:
        return None
    return round(abs(score) / hac, 4)


def le_grade(le: float | None) -> str:
    if le is None:
        return "N/A"
    if le >= 0.4:
        return "◎ 優秀"
    if le >= 0.3:
        return "○ 良好"
    if le >= 0.2:
        return "△ 許容"
    return "✕ 非効率"


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
                    records.append({
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
                    })
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
                        records.append({
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
                        })
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

    # ② スコア上位 TOP_SCORE_N 件
    by_score = sorted(valid, key=lambda r: r["score"])
    top_score = by_score[:TOP_SCORE_N]

    # ③ LE 上位 TOP_LE_N 件 (LE が計算できるもの限定)
    with_le = [r for r in valid if r["le"] is not None]
    by_le   = sorted(with_le, key=lambda r: r["le"], reverse=True)
    top_le  = by_le[:TOP_LE_N]

    # ④ 統合 (重複除去、名前で判断)
    seen  = set()
    final = []
    for r in top_score + top_le:
        if r["name"] not in seen:
            seen.add(r["name"])
            final.append(r)

    # ⑤ 最終スコア順でソート
    final.sort(key=lambda r: r["score"])
    return final


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
        f"| スコア上位 | Top {TOP_SCORE_N} (より負 = 強い結合) |",
        f"| LE 上位 | Top {TOP_LE_N} (重複除去後に追加) |",
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
        f"| 順位 | スコア (kcal/mol) | HAC | LE | LE評価 | 対象残基ペア | 分子名 |",
        f"|------|-----------------|-----|-----|--------|------------|--------|",
    ]

    for rank, r in enumerate(selected, 1):
        le_str = f"{r['le']:.4f}" if r.get("le") is not None else "N/A"
        mw_str = f"{float(r['mw']):.1f}" if r.get("mw") else "?"
        lines.append(
            f"| {rank} | {r['score']:.3f} | {r.get('hac','?')} | {le_str} "
            f"| {le_grade(r.get('le'))} | {r.get('pair','')} | {r['name'][:50]} |"
        )

    lines += [
        "",
        "## ファイル構成",
        "```",
        "Result_Best/",
        "├── 01_score*_<name>.sdf   # ドッキングポーズ SDF (PyMOL/VMD で可視化可)",
        "├── ...                    # 02〜N 番まで同様",
        "├── summary.csv            # 全選抜分子の詳細データ",
        "└── README.md              # 本ファイル",
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
# メイン
# ──────────────────────────────────────────────────────────

def main():
    print("\n" + "=" * 60)
    print("  Result_Best 選抜スクリプト")
    print("=" * 60)

    # 出力ディレクトリ作成
    OUT_DIR.mkdir(exist_ok=True)
    print(f"  出力先: {OUT_DIR.resolve()}\n")

    # 全結果読み込み
    print("  ① 全ドッキング結果を読み込み中...")
    records = load_all_results()
    print(f"     合計 {len(records)} 件 (ドッキングスコアあり)")

    # 選抜
    print(f"\n  ② 選抜 (DrugLike + スコア ≤ {SCORE_CUTOFF} + Top{TOP_SCORE_N} + LE Top{TOP_LE_N})")
    selected = select_best(records)
    print(f"     → {len(selected)} 件 選抜\n")

    if not selected:
        print("  選抜分子なし。基準を緩めるか results/ を確認してください。")
        return

    # SDF コピー
    print("  ③ SDF ファイルをコピー中...")
    selected = copy_sdfs(selected, OUT_DIR)

    # サマリー CSV
    write_summary_csv(selected, OUT_DIR)

    # README
    write_readme(selected, OUT_DIR)

    # ─── 最終サマリー表示 ───
    print()
    print("=" * 90)
    print(f"  {'順位':<5} {'スコア':>9}  {'HAC':>5}  {'LE':>7}  {'LE評価':<10}  分子名")
    print("=" * 90)
    for rank, r in enumerate(selected, 1):
        le_str  = f"{r['le']:.4f}" if r.get("le") is not None else "  N/A "
        grade   = le_grade(r.get("le"))
        src_tag = "STEP3" if "STEP3" in r["source"] else "STEP7"
        print(f"  {rank:<5} {r['score']:>9.3f}  {str(r.get('hac','?')):>5}  {le_str:>7}  "
              f"{grade:<10}  [{src_tag}] {r['name'][:45]}")
    print("=" * 90)
    print(f"\n  完了! Result_Best/ に {len(selected)} 件のSDF + summary.csv + README.md を保存しました")


if __name__ == "__main__":
    main()
