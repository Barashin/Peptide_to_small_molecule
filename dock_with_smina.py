"""
dock_with_smina.py
==================
smina を使ってペプチドの結合部位に低分子候補をドッキングするモジュール。

手順:
  1. タンパク質チェーン (Chain A) を receptor.pdb として抽出
  2. ペプチドチェーン (Chain B) を autobox_ligand 用の参照リガンドとして抽出
  3. smina でドッキング実行 (--autobox_ligand でボックス自動定義)
  4. 結果をスコア順に並べ、CSV・可視化を出力
"""

import os
import re
import subprocess
import json
import csv
import numpy as np
from pathlib import Path


SMINA_PATH = os.path.join(os.path.dirname(__file__), "smina.osx.12")


# ────────────────────────────────────────────
# 受容体 / 参照リガンドの準備
# ────────────────────────────────────────────

def extract_chain(pdb_path: str, chain_id: str, output_path: str):
    """PDB から指定チェーンを抽出して保存"""
    lines = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")):
                if line[21] == chain_id:
                    lines.append(line)
            elif line.startswith("TER"):
                if line[21:22].strip() == chain_id or line[21:22] == chain_id:
                    lines.append(line)
    lines.append("END\n")
    with open(output_path, "w") as f:
        f.writelines(lines)
    print(f"  チェーン {chain_id} を保存: {output_path}  ({len(lines)-1} 行)")


def get_peptide_center(pdb_path: str, chain_id: str) -> tuple[float, float, float]:
    """ペプチドチェーンの重心座標を計算"""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[21] == chain_id:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except ValueError:
                    continue
    if not coords:
        raise ValueError(f"チェーン {chain_id} の座標が見つかりません")
    arr = np.array(coords)
    center = arr.mean(axis=0)
    return float(center[0]), float(center[1]), float(center[2])


def get_peptide_bbox(pdb_path: str, chain_id: str,
                     buffer: float = 6.0) -> dict:
    """ペプチドの座標範囲からドッキングボックスを計算"""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and line[21] == chain_id:
                try:
                    x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                    coords.append((x, y, z))
                except ValueError:
                    continue
    arr = np.array(coords)
    mn, mx = arr.min(axis=0), arr.max(axis=0)
    center = ((mn + mx) / 2).tolist()
    size   = ((mx - mn) + buffer * 2).tolist()
    # 最小20Å確保
    size   = [max(s, 20.0) for s in size]
    return {
        "center_x": round(center[0], 2),
        "center_y": round(center[1], 2),
        "center_z": round(center[2], 2),
        "size_x"  : round(size[0], 1),
        "size_y"  : round(size[1], 1),
        "size_z"  : round(size[2], 1),
    }


# ────────────────────────────────────────────
# smina ドッキング実行
# ────────────────────────────────────────────

def run_smina(receptor_pdb: str,
              ligand_sdf : str,
              peptide_ref: str,
              output_sdf : str,
              log_path   : str,
              exhaustiveness: int = 8,
              num_modes  : int = 9,
              autobox_add: float = 4.0,
              extra_args : list[str] | None = None) -> str:
    """
    smina を実行してドッキングポーズを生成。
    --autobox_ligand でペプチドの結合部位を自動定義。
    Returns: smina の stdout
    """
    cmd = [
        SMINA_PATH,
        "--receptor",      receptor_pdb,
        "--ligand",        ligand_sdf,
        "--autobox_ligand",peptide_ref,
        "--autobox_add",   str(autobox_add),
        "--out",           output_sdf,
        "--log",           log_path,
        "--exhaustiveness",str(exhaustiveness),
        "--num_modes",     str(num_modes),
        "--energy_range",  "5",
    ]
    if extra_args:
        cmd.extend(extra_args)

    print(f"\n  smina 実行中...")
    print(f"  受容体 : {receptor_pdb}")
    print(f"  リガンド: {ligand_sdf}")
    print(f"  参照   : {peptide_ref} (autobox)")
    print(f"  出力   : {output_sdf}")
    print()

    result = subprocess.run(cmd, capture_output=True, text=True)
    stdout = result.stdout + result.stderr

    if result.returncode != 0:
        print(f"  [警告] smina が非ゼロ終了: {result.returncode}")
        print(stdout[-2000:])
    return stdout


def run_smina_per_molecule(receptor_pdb: str,
                           ligand_sdf  : str,
                           peptide_ref : str,
                           output_dir  : str,
                           exhaustiveness: int = 8,
                           num_modes   : int = 9) -> list[dict]:
    """
    SDF 中の各分子を個別にドッキングしてスコアを収集。
    Returns: list of {name, smiles, score, output_sdf}
    """
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(ligand_sdf, removeHs=False)
    results = []

    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"mol_{i}"
        smiles = mol.GetProp("SMILES") if mol.HasProp("SMILES") else ""

        # 個別SDFに書き出し
        single_sdf = os.path.join(output_dir, f"input_{i}_{name[:20]}.sdf")
        out_sdf    = os.path.join(output_dir, f"docked_{i}_{name[:20]}.sdf")
        log_file   = os.path.join(output_dir, f"log_{i}_{name[:20]}.txt")

        w = Chem.SDWriter(single_sdf)
        w.write(mol)
        w.close()

        stdout = run_smina(receptor_pdb, single_sdf, peptide_ref,
                           out_sdf, log_file,
                           exhaustiveness=exhaustiveness,
                           num_modes=num_modes)

        # ログからスコアを解析
        scores = parse_smina_scores(log_file, stdout)
        best   = min(scores) if scores else None

        results.append({
            "index"     : i,
            "name"      : name,
            "smiles"    : smiles,
            "best_score": best,
            "all_scores": scores,
            "output_sdf": out_sdf if os.path.exists(out_sdf) else None,
            "log"       : log_file,
        })

        flag = f"{best:.3f} kcal/mol" if best is not None else "スコアなし"
        print(f"  [{i+1}] {name:<35}  Best: {flag}")

    return results


# ────────────────────────────────────────────
# ログ解析
# ────────────────────────────────────────────

def parse_smina_scores(log_path: str, stdout: str = "") -> list[float]:
    """smina ログからアフィニティスコアを抽出"""
    scores = []
    # ログファイルから
    text = ""
    if os.path.exists(log_path):
        with open(log_path) as f:
            text = f.read()
    text += "\n" + stdout

    # smina 出力フォーマット:
    #    1       -7.123     0.000     0.000
    pattern = re.compile(r"^\s*\d+\s+([-\d.]+)\s+[\d.]+\s+[\d.]+", re.MULTILINE)
    for m in pattern.finditer(text):
        try:
            scores.append(float(m.group(1)))
        except ValueError:
            pass
    return scores


# ────────────────────────────────────────────
# smina --score_only リスコアリング
# ────────────────────────────────────────────

def parse_score_only_output(stdout: str) -> dict | None:
    """
    smina --score_only の出力からスコア項を解析する。

    smina 実際の出力フォーマット:
      Weights      Terms
      -0.035579    gauss(o=0,_w=0.5,_c=8)
      -0.005156    gauss(o=3,_w=2,_c=8)
       0.840245    repulsion(o=0,_c=8)
      -0.035069    hydrophobic(g=0.5,_b=1.5,_c=8)
      -0.587439    non_dir_h_bond(g=-0.7,_b=0,_c=8)

      ## Name gauss(o=0,...) gauss(o=3,...) repulsion(...) hydrophobic(...) non_dir_h_bond(...) num_tors_div
      Affinity: -5.514 (kcal/mol)
      Intramolecular energy: 0.000
      Term values, before weighting:
      ## mol_name  54.686  534.676  0.908  44.905  0.000  0.000

    スコア項はraw値 × 重みで重みつき寄与を計算する。
    """
    result = {}

    # アフィニティ
    affinity_m = re.search(r"Affinity:\s*([-\d.]+)", stdout)
    if affinity_m:
        result["affinity"] = float(affinity_m.group(1))

    intra_m = re.search(r"Intramolecular energy:\s*([-\d.]+)", stdout)
    if intra_m:
        result["intramolecular"] = float(intra_m.group(1))

    # 重みと項名を抽出
    # "Weights      Terms" セクション: 各行 "  重み  項名"
    weights_map = {}
    weight_pattern = re.compile(
        r"^\s*([-\d.]+)\s+(gauss\([^)]+\)|repulsion\([^)]+\)|"
        r"hydrophobic\([^)]+\)|non_dir_h_bond\([^)]+\))",
        re.MULTILINE
    )
    for m in weight_pattern.finditer(stdout):
        weights_map[m.group(2)] = float(m.group(1))

    # ヘッダ行: ## Name term1 term2 ...
    header_m = re.search(r"^## Name\s+(.+)$", stdout, re.MULTILINE)
    # raw値行: ## mol_name val1 val2 ...
    values_m = re.search(r"^## (?!Name)(\S+)\s+([\d.\s-]+)$", stdout, re.MULTILINE)

    if header_m and values_m:
        headers = header_m.group(1).split()
        raw_vals = [float(v) for v in values_m.group(2).split()]

        # 各項の重みつき寄与を計算
        weighted = {}
        for i, (hdr, raw) in enumerate(zip(headers, raw_vals)):
            # ヘッダ名に一致する重みを探す
            w = None
            for term_name, term_w in weights_map.items():
                # ヘッダはsmall shortname, term_nameはfull name で照合
                if ("gauss" in hdr and "gauss" in term_name):
                    # gauss(o=0...) → gauss1, gauss(o=3...) → gauss2 で区別
                    if i == 0 and "o=0" in term_name:
                        w = term_w
                    elif i == 1 and "o=3" in term_name:
                        w = term_w
                elif "repulsion" in hdr and "repulsion" in term_name:
                    w = term_w
                elif "hydrophobic" in hdr and "hydrophobic" in term_name:
                    w = term_w
                elif "h_bond" in hdr and "h_bond" in term_name:
                    w = term_w
            if w is not None:
                weighted[hdr] = round(w * raw, 5)

        # 統一キー名でマップ
        for hdr, val in weighted.items():
            if "gauss" in hdr and "gauss1" not in result:
                result["gauss1"] = val
            elif "gauss" in hdr:
                result["gauss2"] = val
            elif "repulsion" in hdr:
                result["repulsion"] = val
            elif "hydrophobic" in hdr:
                result["hydrophobic"] = val
            elif "h_bond" in hdr:
                result["Hbond"] = val

    return result if result else None


def rescore_with_smina(receptor_pdb: str,
                       docked_sdf  : str,
                       output_dir  : str,
                       mol_name    : str = "mol") -> list[dict]:
    """
    ドッキングポーズSDFを --score_only で再スコアリングし、
    Vinaスコア項の分解（疎水性・H結合・立体反発など）を取得する。

    これにより、どの相互作用が結合に貢献しているかを定量的に把握できる。
    （MM-GBSAのエネルギー分解に相当する簡易版）

    Args:
        receptor_pdb: 受容体PDBパス
        docked_sdf  : ドッキングポーズSDF（複数ポーズ含む）
        output_dir  : 出力ディレクトリ
        mol_name    : 分子名（ファイル名用）
    Returns:
        list of {pose_id, affinity, gauss1, gauss2, repulsion, hydrophobic, Hbond}
    """
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(docked_sdf, removeHs=False)
    pose_scores = []

    for pose_id, mol in enumerate(suppl):
        if mol is None:
            continue

        # ポーズを個別SDFに書き出し
        pose_sdf = os.path.join(output_dir,
                                f"pose_{mol_name[:15]}_{pose_id}.sdf")
        w = Chem.SDWriter(pose_sdf)
        w.write(mol)
        w.close()

        cmd = [
            SMINA_PATH,
            "--receptor",   receptor_pdb,
            "--ligand",     pose_sdf,
            "--score_only",
            "--log",        os.path.join(output_dir,
                                         f"rescore_{mol_name[:15]}_{pose_id}.txt"),
        ]
        result = subprocess.run(cmd, capture_output=True, text=True)
        stdout = result.stdout + result.stderr
        parsed = parse_score_only_output(stdout)

        if parsed:
            parsed["pose_id"] = pose_id
            pose_scores.append(parsed)

    return pose_scores


def rescore_all_docked(docking_results: list[dict],
                       receptor_pdb   : str,
                       output_dir     : str) -> list[dict]:
    """
    run_smina_per_molecule() の結果に対して全ポーズをリスコアリングし、
    元の結果dictにスコア分解を追加して返す。
    """
    rescore_dir = os.path.join(output_dir, "rescoring")
    os.makedirs(rescore_dir, exist_ok=True)

    enriched = []
    for r in docking_results:
        if r["output_sdf"] is None or not os.path.exists(r["output_sdf"]):
            enriched.append({**r, "score_terms": []})
            continue

        terms = rescore_with_smina(
            receptor_pdb, r["output_sdf"], rescore_dir,
            mol_name=r["name"][:15]
        )
        enriched.append({**r, "score_terms": terms})

        if terms:
            best = min(terms, key=lambda x: x.get("affinity", 0))
            print(f"  [{r['name'][:30]}] rescore best: "
                  f"{best.get('affinity', 'N/A'):.3f} kcal/mol  "
                  f"(HB={best.get('Hbond', 0):.2f}, "
                  f"HYD={best.get('hydrophobic', 0):.2f}, "
                  f"REP={best.get('repulsion', 0):.2f})")

    return enriched


def plot_score_decomposition(enriched_results: list[dict], output_path: str):
    """
    スコア項分解のスタック棒グラフ。
    各分子のベストポーズについて gauss / hydrophobic / Hbond / repulsion を可視化。
    """
    import plot_utils  # noqa: F401
    import matplotlib.pyplot as plt
    import numpy as np

    records = []
    for r in enriched_results:
        if not r.get("score_terms"):
            continue
        best = min(r["score_terms"], key=lambda x: x.get("affinity", 0))
        records.append({
            "name"       : r["name"][:22],
            "gauss"      : best.get("gauss1", 0) + best.get("gauss2", 0),
            "hydrophobic": best.get("hydrophobic", 0),
            "hbond"      : best.get("Hbond", 0),
            "repulsion"  : best.get("repulsion", 0),
            "affinity"   : best.get("affinity", 0),
        })

    if not records:
        return

    names = [r["name"] for r in records]
    x     = np.arange(len(names))
    width = 0.6

    fig, ax = plt.subplots(figsize=(max(8, len(names) * 1.5), 6))

    # 負のスコア（有利）を下方向に、正（不利）を上方向に積み上げ
    neg_terms = [("gauss",       "#3498db", "Gauss (vdW)"),
                 ("hydrophobic", "#f1c40f", "疎水性"),
                 ("hbond",       "#e74c3c", "H結合")]
    pos_terms = [("repulsion",   "#95a5a6", "立体反発")]

    bottom_neg = np.zeros(len(names))
    for key, color, label in neg_terms:
        vals = np.array([r[key] for r in records])
        ax.bar(x, vals, width, bottom=bottom_neg,
               color=color, label=label, edgecolor="white")
        bottom_neg += vals

    bottom_pos = np.zeros(len(names))
    for key, color, label in pos_terms:
        vals = np.array([r[key] for r in records])
        ax.bar(x, vals, width, bottom=bottom_pos,
               color=color, label=label, edgecolor="white", alpha=0.7)
        bottom_pos += vals

    # 合計アフィニティをプロット
    affinities = [r["affinity"] for r in records]
    ax.plot(x, affinities, "ko--", markersize=7, linewidth=1.5,
            label="合計 Affinity", zorder=5)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.set_ylabel("スコア寄与 (kcal/mol)", fontsize=11)
    ax.set_title("smina スコア項分解 (ベストポーズ)", fontsize=13, fontweight="bold")
    ax.legend(loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  スコア分解グラフ 保存: {output_path}")


def save_rescore_csv(enriched_results: list[dict], output_path: str):
    import csv
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["name", "pose_id", "affinity",
                         "gauss1", "gauss2", "repulsion",
                         "hydrophobic", "Hbond"])
        for r in enriched_results:
            for t in r.get("score_terms", []):
                writer.writerow([
                    r["name"], t.get("pose_id", ""),
                    t.get("affinity", ""), t.get("gauss1", ""),
                    t.get("gauss2", ""),  t.get("repulsion", ""),
                    t.get("hydrophobic", ""), t.get("Hbond", ""),
                ])
    print(f"  rescoring CSV 保存: {output_path}")


# ────────────────────────────────────────────
# 結果の保存・表示
# ────────────────────────────────────────────

def save_docking_results(results: list[dict], output_path: str):
    """ドッキング結果をCSVに保存"""
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["rank", "name", "best_score_kcal/mol",
                         "all_scores", "smiles", "output_sdf"])
        ranked = sorted(
            [r for r in results if r["best_score"] is not None],
            key=lambda x: x["best_score"]
        )
        # スコアなし分子を末尾に追加
        ranked += [r for r in results if r["best_score"] is None]

        for rank, r in enumerate(ranked, 1):
            writer.writerow([
                rank,
                r["name"],
                r["best_score"] if r["best_score"] is not None else "N/A",
                ";".join(f"{s:.3f}" for s in r["all_scores"]),
                r["smiles"],
                r["output_sdf"] or "",
            ])
    print(f"  ドッキング結果 CSV 保存: {output_path}")


def print_docking_summary(results: list[dict]):
    """ドッキング結果をターミナルに表示"""
    ranked = sorted(
        [r for r in results if r["best_score"] is not None],
        key=lambda x: x["best_score"]
    )
    print("\n" + "=" * 65)
    print("  ドッキング結果 ランキング")
    print("=" * 65)
    print(f"  {'順位':<5} {'分子名':<35} {'Best Score (kcal/mol)'}")
    print("  " + "-" * 60)
    for rank, r in enumerate(ranked, 1):
        bar = "█" * int(abs(r["best_score"]) * 2)
        print(f"  {rank:<5} {r['name']:<35} {r['best_score']:>8.3f}  {bar}")
    no_score = [r for r in results if r["best_score"] is None]
    for r in no_score:
        print(f"  {'--':<5} {r['name']:<35} {'N/A':>8}")
    print("=" * 65)


def plot_docking_scores(results: list[dict], output_path: str):
    """ドッキングスコアの棒グラフ"""
    import plot_utils  # noqa: F401
    import matplotlib.pyplot as plt

    ranked = sorted(
        [r for r in results if r["best_score"] is not None],
        key=lambda x: x["best_score"]
    )
    if not ranked:
        return

    names  = [r["name"][:25] for r in ranked]
    scores = [r["best_score"] for r in ranked]
    colors = ["#e74c3c" if s < -7 else "#e67e22" if s < -5 else "#3498db"
              for s in scores]

    fig, ax = plt.subplots(figsize=(max(8, len(names) * 1.4), 5))
    bars = ax.barh(names[::-1], [abs(s) for s in scores[::-1]],
                   color=colors[::-1], edgecolor="white")
    ax.set_xlabel("|Affinity| (kcal/mol) — 右ほど強い結合", fontsize=11)
    ax.set_title("smina ドッキングスコア比較", fontsize=13, fontweight="bold")

    for bar, s in zip(bars, scores[::-1]):
        ax.text(bar.get_width() + 0.05, bar.get_y() + bar.get_height() / 2,
                f"{s:.2f}", va="center", fontsize=9)

    import matplotlib.patches as mpatches
    patches = [
        mpatches.Patch(color="#e74c3c", label="< -7 (強い)"),
        mpatches.Patch(color="#e67e22", label="-7 ~ -5 (中程度)"),
        mpatches.Patch(color="#3498db", label="> -5 (弱い)"),
    ]
    ax.legend(handles=patches, loc="lower right", fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  ドッキングスコア図 保存: {output_path}")


# ────────────────────────────────────────────
# メイン (単体実行)
# ────────────────────────────────────────────

if __name__ == "__main__":
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("pdb",        nargs="?", default="Protein_Peptide.pdb")
    p.add_argument("--ligands",  default="results/candidate_ligands.sdf")
    p.add_argument("--output-dir", default="results/docking")
    p.add_argument("--exhaustiveness", type=int, default=8)
    p.add_argument("--num-modes",      type=int, default=9)
    args = p.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # チェーン抽出
    receptor_pdb = os.path.join(args.output_dir, "receptor.pdb")
    peptide_ref  = os.path.join(args.output_dir, "peptide_ref.pdb")
    extract_chain(args.pdb, "A", receptor_pdb)
    extract_chain(args.pdb, "B", peptide_ref)

    # ドッキングボックス確認
    box = get_peptide_bbox(args.pdb, "B")
    print(f"\n  ドッキングボックス: center({box['center_x']}, {box['center_y']}, "
          f"{box['center_z']})  size({box['size_x']} x {box['size_y']} x {box['size_z']})")

    # ドッキング実行
    results = run_smina_per_molecule(
        receptor_pdb, args.ligands, peptide_ref,
        args.output_dir,
        exhaustiveness=args.exhaustiveness,
        num_modes=args.num_modes,
    )

    print_docking_summary(results)
    save_docking_results(results, os.path.join(args.output_dir, "docking_results.csv"))
    plot_docking_scores(results, os.path.join(args.output_dir, "docking_scores.png"))

    with open(os.path.join(args.output_dir, "docking_results.json"), "w") as f:
        json.dump(results, f, indent=2, ensure_ascii=False)
