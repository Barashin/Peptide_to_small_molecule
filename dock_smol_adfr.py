#!/usr/bin/env python3
"""
dock_smol_adfr.py
=================
AutoDockFR (ADFR) を使って低分子候補をドッキングし、
AutoDock CrankPep (ADCP) による環状ペプチド結果と比較する。

【ADFR の位置づけ】
  ADFR と ADCP は同じ AutoDock 力場ベース:
    ✓ 静電相互作用 (クーロン則)
    ✓ 脱溶媒和項 (原子ベース)
    ✓ vdW / HBond
    ✓ 同じ .trg ターゲットファイルを共用
  smina(Vina) は経験的スコア (静電なし、溶媒なし)

【ADFR と ADCP の残る違い】
  ADCP: AutoDock 力場 + GB/SA implicit 溶媒 → 典型 -10〜-35 kcal/mol
  ADFR: AutoDock 力場のみ (implicit 溶媒なし) → 典型  -5〜-15 kcal/mol
  → smina よりは比較しやすいが、数値の直接比較は引き続き注意が必要

【前提条件】
  1. adcpsuite 環境のセットアップ:
       bash adcpsuite_micromamba.sh
  2. receptor.trg の作成 (ADCP と共用):
       python dock_cyclic_adcp.py --setup-receptor
  3. 低分子候補の選抜 (任意、推奨):
       python collect_best.py

Usage:
  # 通常実行 (Result_Best の上位 15 件を ADFR でドッキング)
  python dock_smol_adfr.py

  # クイックテスト (N=5, n=500000, ≈1〜2分/分子)
  python dock_smol_adfr.py --quick

  # 処理件数を絞る
  python dock_smol_adfr.py --top-n 5

  # 既存結果を上書きして再実行
  python dock_smol_adfr.py --force

出力:
  results/adfr_docking/
  ├── ligands_pdbqt/               # SDF → PDBQT 変換済みファイル
  ├── {mol_name}/                  # 各分子の ADFR 出力
  │   ├── result_{name}_summary.dlg
  │   └── result_{name}_out.pdbqt
  ├── adfr_docking_scores.png      # ADFR スコア棒グラフ (LE 色分け)
  ├── adfr_vs_smina.png            # 同一分子: ADFR vs smina 並列比較
  ├── adfr_vs_adcp_comparison.png  # ADFR 低分子 vs ADCP 環状ペプチド
  ├── le_comparison_3methods.png   # LE: ADCP / ADFR / smina 3系統比較
  ├── adfr_docking_results.csv     # 全結果 CSV
  └── adfr_comparison_report.json  # 結果サマリー JSON
"""

import os
import re
import csv
import sys
import json
import shutil
import subprocess
import argparse
from pathlib import Path
from dataclasses import dataclass, field
from typing import Optional, List, Tuple

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


# ──────────────────────────────────────────────────────────────────────────────
# 定数
# ──────────────────────────────────────────────────────────────────────────────

BASE_DIR    = Path(__file__).parent
RESULTS_DIR = BASE_DIR / "results"
ADCP_DIR    = RESULTS_DIR / "adcp_docking"
ADFR_DIR    = RESULTS_DIR / "adfr_docking"

DEFAULT_N_RUNS  = 20
DEFAULT_N_EVALS = 2_500_000
QUICK_N_RUNS    = 5
QUICK_N_EVALS   = 500_000


# ──────────────────────────────────────────────────────────────────────────────
# データクラス
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class ADFRResult:
    """1 分子の ADFR ドッキング結果"""
    name:          str
    affinity:      float          # Best affinity mode 1 (kcal/mol)
    best_energy:   float          # Best MC/GA energy
    n_clusters:    int
    status:        str            # success / failed / timeout / convert_failed / error
    error_message: Optional[str] = None
    dlg_file:      Optional[str] = None
    pose_file:     Optional[str] = None   # .pdbqt or .pdb
    hac:           Optional[int] = None
    le_adfr:       Optional[float] = None
    smina_score:   Optional[float] = None  # 同一分子の smina スコア (3-way 比較用)
    smina_le:      Optional[float] = None
    pair:          str = ""


# ──────────────────────────────────────────────────────────────────────────────
# 環境チェック (dock_cyclic_adcp.py と同様)
# ──────────────────────────────────────────────────────────────────────────────

def find_micromamba() -> Optional[str]:
    """micromamba の実行パスを探す"""
    try:
        r = subprocess.run(["which", "micromamba"], capture_output=True, text=True)
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass
    for p in [
        os.path.expanduser("~/.local/bin/micromamba"),
        os.path.expanduser("~/micromamba/bin/micromamba"),
        os.path.expanduser("~/mambaforge/bin/micromamba"),
        "/usr/local/bin/micromamba",
        "/opt/homebrew/bin/micromamba",
    ]:
        if os.path.exists(p):
            return p
    return None


def find_adcpsuite_env(mamba: str) -> Optional[str]:
    """adcpsuite 環境のフルパスを返す（異なるconda/mambaルートにも対応）"""
    try:
        r = subprocess.run([mamba, "env", "list"], capture_output=True, text=True)
        for line in r.stdout.splitlines():
            if "adcpsuite" in line:
                parts = line.split()
                path = parts[-1]
                if os.path.isdir(path):
                    return path
    except Exception:
        pass
    for p in [
        os.path.expanduser("~/miniforge3/envs/adcpsuite"),
        os.path.expanduser("~/miniforge/envs/adcpsuite"),
        os.path.expanduser("~/mambaforge/envs/adcpsuite"),
        os.path.expanduser("~/micromamba/envs/adcpsuite"),
        os.path.expanduser("~/miniconda3/envs/adcpsuite"),
        os.path.expanduser("~/opt/anaconda3/envs/adcpsuite"),
        "/opt/anaconda3/envs/adcpsuite",
        "/opt/homebrew/Caskroom/miniforge/base/envs/adcpsuite",
    ]:
        if os.path.isdir(p):
            return p
    return None


def check_adcpsuite_env(mamba: str) -> bool:
    return find_adcpsuite_env(mamba) is not None


def run_in_adcpsuite(mamba: str, cmd: str, cwd: str = None,
                      timeout: int = 600) -> subprocess.CompletedProcess:
    """adcpsuite 環境でシェルコマンドを実行（-p フルパス優先）"""
    env_path = find_adcpsuite_env(mamba)
    if env_path:
        full_cmd = f"{mamba} run -p {env_path} {cmd}"
    else:
        full_cmd = f"{mamba} run -n adcpsuite {cmd}"
    return subprocess.run(
        full_cmd, shell=True, capture_output=True, text=True,
        cwd=cwd, timeout=timeout,
    )


# ──────────────────────────────────────────────────────────────────────────────
# ターゲットファイル確認
# ──────────────────────────────────────────────────────────────────────────────

def find_receptor_trg() -> Optional[Path]:
    """dock_cyclic_adcp.py が作成した receptor.trg を探す"""
    trg = ADCP_DIR / "receptor.trg"
    return trg if trg.exists() else None


# ──────────────────────────────────────────────────────────────────────────────
# SDF → PDBQT 変換
# ──────────────────────────────────────────────────────────────────────────────

def convert_sdf_to_pdbqt(mamba: str, sdf_path: Path,
                          out_pdbqt: Path, force: bool = False) -> bool:
    """
    SDF ファイルを AutoDock PDBQT 形式に変換する。

    obabel (adcpsuite 環境) を使用し、Gasteiger 電荷を付与する。
    obabel が失敗した場合は prepare_ligand4.py (ADFRsuite) にフォールバック。
    """
    if out_pdbqt.exists() and not force:
        return True

    out_pdbqt.parent.mkdir(parents=True, exist_ok=True)

    # obabel で変換 (水素付加 + Gasteiger 電荷)
    cmd = (
        f'obabel -isdf "{sdf_path.resolve()}"'
        f' -opdbqt -O "{out_pdbqt.resolve()}"'
        f' --partialcharge gasteiger -h'
    )
    r = run_in_adcpsuite(mamba, cmd, timeout=60)

    if r.returncode != 0 or not out_pdbqt.exists():
        # フォールバック: -h なしで試みる
        cmd = (
            f'obabel -isdf "{sdf_path.resolve()}"'
            f' -opdbqt -O "{out_pdbqt.resolve()}"'
            f' --partialcharge gasteiger'
        )
        r = run_in_adcpsuite(mamba, cmd, timeout=60)

    if r.returncode != 0 or not out_pdbqt.exists():
        # 最終フォールバック: prepare_ligand4.py
        cmd = (
            f'prepare_ligand4.py'
            f' -l "{sdf_path.resolve()}"'
            f' -o "{out_pdbqt.resolve()}"'
        )
        r = run_in_adcpsuite(mamba, cmd, timeout=60)

    # obabel は複数コンフォーマー(MODEL)を出力するが ADFR は単一モデルを期待。
    # 最初の MODEL だけ抽出し、ENDMDL / MODEL 行を除去する。
    if out_pdbqt.exists():
        lines = out_pdbqt.read_text().splitlines()
        first_model = []
        in_first = True
        model_count = 0
        for l in lines:
            if l.strip().startswith("MODEL"):
                model_count += 1
                if model_count > 1:
                    in_first = False
                continue  # MODEL 行自体は除去
            if l.strip() == "ENDMDL":
                if in_first:
                    continue  # 除去
                break  # 2番目以降のモデルは無視
            if in_first:
                first_model.append(l)
        out_pdbqt.write_text("\n".join(first_model) + "\n")

    return out_pdbqt.exists()


# ──────────────────────────────────────────────────────────────────────────────
# DLG ファイルのパース (ADFR / ADCP 共通フォーマット)
# ──────────────────────────────────────────────────────────────────────────────

def parse_dlg_file(dlg_file: Path) -> Tuple[float, float, int]:
    """
    AutoDock .dlg ファイルから affinity, best_energy, n_clusters を抽出。

    ADFR 形式 (NoName*.dlg):
        _Gen0044 Score: -44.972 LL: -6.677 LR: -38.295 evals: 255101  2
        CNUM  len best  Rmsd   Score    FEB    <Score>  stdev cluster
          0    1    0  -1.00  -44.972  -38.295  -44.972  0.000 [0]
        FEB = Free Energy of Binding (≈ affinity)
        Score = 内部スコア (≈ best_energy)

    ADCP 形式 (旧フォーマットも対応):
        mode |  affinity  | ...
           1      -8.50   ...
        bestEnergy in run 3 -9.12 (0)
    """
    affinity    = float("inf")
    best_energy = float("inf")
    n_clusters  = 0

    if not dlg_file.exists():
        return affinity, best_energy, n_clusters

    content = dlg_file.read_text(errors="replace")

    # --- ADFR 形式: _GenNNNN Score: ... LR: <FEB> ---
    gen_matches = re.findall(
        r"_Gen\d+\s+Score:\s+([-\d\.]+)\s+LL:\s+[-\d\.]+\s+LR:\s+([-\d\.]+)",
        content)
    if gen_matches:
        # 最後の _Gen 行が最良結果
        best_score_str, best_feb_str = gen_matches[-1]
        affinity    = float(best_feb_str)   # FEB = LR
        best_energy = float(best_score_str) # Score

        # cluster 数: CNUM テーブルの行数 (最後のクラスタリングブロック)
        cluster_lines = re.findall(
            r"^\s+\d+\s+\d+\s+\d+\s+[-\d\.]+\s+[-\d\.]+\s+[-\d\.]+",
            content, re.MULTILINE)
        n_clusters = len(cluster_lines) if cluster_lines else 1
        return affinity, best_energy, n_clusters

    # --- ADCP 形式 (旧フォーマット) ---
    m = re.search(r"^\s+1\s+([-\d\.]+)", content, re.MULTILINE)
    if m:
        affinity = float(m.group(1))
    m = re.search(r"bestEnergy in run \d+ ([-\d\.]+)", content)
    if m:
        best_energy = float(m.group(1))
    n_clusters = len(re.findall(r"^\s+\d+\s+[-\d\.]+", content, re.MULTILINE))

    return affinity, best_energy, n_clusters


def _find_dlg(out_dir: Path, job_name: str) -> Optional[Path]:
    """ADFR の出力 DLG ファイルを探す。

    ADFR は result_<job>/NoName*.dlg に個別結果を出力するため、
    その中で最もサイズの大きいもの (最終 GA run) を返す。
    見つからない場合は summary.dlg にフォールバック。
    """
    # ADFR: サブディレクトリ内の NoName*.dlg
    sub_dir = out_dir / f"result_{job_name}"
    if sub_dir.is_dir():
        dlg_files = sorted(sub_dir.glob("*.dlg"),
                           key=lambda p: p.stat().st_size, reverse=True)
        if dlg_files:
            return dlg_files[0]

    # フォールバック: summary.dlg
    expected = out_dir / f"result_{job_name}_summary.dlg"
    if expected.exists():
        return expected
    dlg_files = list(out_dir.glob("*.dlg"))
    return dlg_files[0] if dlg_files else None


def _find_pose(out_dir: Path, job_name: str) -> Optional[Path]:
    """最良ポーズファイル (.pdbqt または .pdb) を探す"""
    for pattern in [f"result_{job_name}_out.pdbqt",
                    f"result_{job_name}_out.pdb",
                    f"*_out.pdbqt", f"*_out.pdb"]:
        matches = list(out_dir.glob(pattern))
        if matches:
            return matches[0]
    return None


# ──────────────────────────────────────────────────────────────────────────────
# ADFR ドッキング (単一分子)
# ──────────────────────────────────────────────────────────────────────────────

def run_adfr_docking(
    mamba:       str,
    name:        str,
    pdbqt_path:  Path,
    trg_file:    Path,
    n_runs:      int,
    n_evals:     int,
    timeout:     int,
    force:       bool,
    hac:         Optional[int]   = None,
    smina_score: Optional[float] = None,
    smina_le:    Optional[float] = None,
    pair:        str             = "",
) -> ADFRResult:
    """
    ADFR (AutoDockFR) で低分子をドッキングする。

    Parameters
    ----------
    name        : 分子名
    pdbqt_path  : 入力 PDBQT ファイル
    trg_file    : 受容体ターゲットファイル (.trg)
    n_runs      : GA 探索回数 (デフォルト: 20)
    n_evals     : 評価ステップ数 (デフォルト: 2,500,000)
    timeout     : タイムアウト秒数
    force       : 既存結果を上書きするか

    ADFRコマンド:
        adfr -t receptor.trg -l ligand.pdbqt
             -o result_name -N 20 -n 2500000
             -w output_dir/
    """
    # ジョブ名: 英数字とハイフン/アンダースコアのみ
    job_name = re.sub(r"[^\w\-]", "_", name)[:50].lower()
    out_dir  = ADFR_DIR / job_name
    out_dir.mkdir(parents=True, exist_ok=True)

    # キャッシュ確認
    dlg_file = _find_dlg(out_dir, job_name)
    if dlg_file and not force:
        print(f"    既存結果を使用: {dlg_file.name}")
        affinity, best_energy, n_clusters = parse_dlg_file(dlg_file)
        pose_file = _find_pose(out_dir, job_name)
        le = round(abs(affinity) / hac, 4) if (hac and affinity < 0) else None
        return ADFRResult(
            name=name, affinity=affinity, best_energy=best_energy,
            n_clusters=n_clusters, status="success",
            dlg_file=str(dlg_file),
            pose_file=str(pose_file) if pose_file else None,
            hac=hac, le_adfr=le,
            smina_score=smina_score, smina_le=smina_le, pair=pair,
        )

    # ADFR コマンド
    # adfr -t receptor.trg -l ligand.pdbqt -o result_job -n 20 -e 2500000
    # 注: -n=nbRuns, -e=maxEvals, -o=logFilename
    cmd = (
        f'adfr'
        f' -t "{trg_file.resolve()}"'
        f' -l "{pdbqt_path.resolve()}"'
        f' -o result_{job_name}'
        f' -n {n_runs}'
        f' -e {n_evals}'
    )

    try:
        r = run_in_adcpsuite(mamba, cmd, cwd=str(out_dir.resolve()),
                             timeout=timeout)

        if r.stderr:
            (out_dir / "adfr_stderr.log").write_text(r.stderr)

        # DLG ファイルが生成されていれば成功とみなす
        # (numpy DeprecationWarning 等で returncode != 0 でも結果は正常な場合がある)
        dlg_file  = _find_dlg(out_dir, job_name)

        if r.returncode != 0 and dlg_file is None:
            return ADFRResult(
                name=name, affinity=float("inf"), best_energy=float("inf"),
                n_clusters=0, status="failed",
                error_message=r.stderr[:400],
                hac=hac, smina_score=smina_score, smina_le=smina_le, pair=pair,
            )
        pose_file = _find_pose(out_dir, job_name)

        if dlg_file is None:
            return ADFRResult(
                name=name, affinity=float("inf"), best_energy=float("inf"),
                n_clusters=0, status="no_output",
                error_message="DLG ファイルが生成されませんでした",
                hac=hac, smina_score=smina_score, smina_le=smina_le, pair=pair,
            )

        affinity, best_energy, n_clusters = parse_dlg_file(dlg_file)
        le = round(abs(affinity) / hac, 4) if (hac and affinity < 0) else None

        return ADFRResult(
            name=name, affinity=affinity, best_energy=best_energy,
            n_clusters=n_clusters, status="success",
            dlg_file=str(dlg_file),
            pose_file=str(pose_file) if pose_file else None,
            hac=hac, le_adfr=le,
            smina_score=smina_score, smina_le=smina_le, pair=pair,
        )

    except subprocess.TimeoutExpired:
        return ADFRResult(
            name=name, affinity=float("inf"), best_energy=float("inf"),
            n_clusters=0, status="timeout",
            error_message=f"Timeout after {timeout}s",
            hac=hac, smina_score=smina_score, smina_le=smina_le, pair=pair,
        )
    except Exception as e:
        return ADFRResult(
            name=name, affinity=float("inf"), best_energy=float("inf"),
            n_clusters=0, status="error", error_message=str(e),
            hac=hac, smina_score=smina_score, smina_le=smina_le, pair=pair,
        )


# ──────────────────────────────────────────────────────────────────────────────
# 分子読み込み
# ──────────────────────────────────────────────────────────────────────────────

def load_molecules(top_n: int = 15) -> List[dict]:
    """
    Result_Best/summary.csv からドッキング対象分子を読み込む。

    SDF ファイルが実在する分子のみを対象とする。
    smina スコアが良い順 (昇順) にソートして返す。
    """
    summary_csv = BASE_DIR / "Result_Best" / "summary.csv"
    if not summary_csv.exists():
        raise FileNotFoundError(
            f"Result_Best/summary.csv が見つかりません。\n"
            f"先に collect_best.py を実行してください:\n"
            f"  python collect_best.py"
        )

    molecules = []
    with open(summary_csv) as f:
        for row in csv.DictReader(f):
            sdf_path = row.get("out_sdf", "").strip()
            if not sdf_path or not Path(sdf_path).exists():
                continue
            try:
                molecules.append({
                    "name":        row["name"],
                    "smina_score": float(row["score"]),
                    "hac":         int(row["hac"])   if row.get("hac")  else None,
                    "smina_le":    float(row["le"])   if row.get("le")   else None,
                    "sdf_path":    Path(sdf_path),
                    "pair":        row.get("pair", ""),
                })
            except (ValueError, KeyError):
                pass

    molecules.sort(key=lambda m: m["smina_score"])
    return molecules[:top_n]


# ──────────────────────────────────────────────────────────────────────────────
# ADCP 結果読み込み
# ──────────────────────────────────────────────────────────────────────────────

def load_adcp_result() -> Optional[dict]:
    """
    dock_cyclic_adcp.py が生成した adcp_comparison_report.json を読み込む。
    ファイルが存在しない場合は None を返す (ADCP 実行前でも動作するよう配慮)。
    """
    report_path = ADCP_DIR / "adcp_comparison_report.json"
    if not report_path.exists():
        return None
    with open(report_path) as f:
        report = json.load(f)
    return report.get("adcp_cyclic_peptide")


# ──────────────────────────────────────────────────────────────────────────────
# 可視化ユーティリティ
# ──────────────────────────────────────────────────────────────────────────────

def _le_color(le: Optional[float]) -> str:
    if le is None:  return "gray"
    if le >= 0.4:   return "#2ecc71"   # 緑: 優秀
    if le >= 0.3:   return "#3498db"   # 青: 良好
    if le >= 0.2:   return "#f39c12"   # 橙: 許容
    return "#e74c3c"                   # 赤: 非効率


def _setup_font():
    try:
        from plot_utils import setup_japanese_font
        setup_japanese_font()
    except ImportError:
        pass


# ──────────────────────────────────────────────────────────────────────────────
# グラフ 1: ADFR スコア棒グラフ
# ──────────────────────────────────────────────────────────────────────────────

def plot_adfr_scores(results: List[ADFRResult], out_dir: Path):
    """ADFR ドッキングスコアを LE 色分き横棒グラフで表示"""
    if not HAS_MATPLOTLIB:
        return
    _setup_font()

    valid = sorted(
        [r for r in results if r.status == "success" and r.affinity < 0],
        key=lambda r: r.affinity,
    )
    if not valid:
        print("  [情報] 成功した ADFR 結果がありません")
        return

    names  = [r.name[:42] + "…" if len(r.name) > 42 else r.name for r in valid]
    scores = [r.affinity for r in valid]
    colors = [_le_color(r.le_adfr) for r in valid]

    fig, ax = plt.subplots(figsize=(14, max(6, len(valid) * 0.5 + 2)))
    bars = ax.barh(range(len(valid)), scores, color=colors, alpha=0.85)
    ax.set_yticks(range(len(valid)))
    ax.set_yticklabels(names, fontsize=8)
    ax.set_xlabel("ADFR Affinity (kcal/mol)  [AutoDock 力場 / 静電+脱溶媒和あり]", fontsize=11)
    ax.set_title(
        "AutoDockFR (ADFR) 低分子ドッキング結果\n"
        "[AutoDock 力場 / implicit 溶媒なし]",
        fontsize=12,
    )
    ax.invert_yaxis()
    ax.invert_xaxis()
    ax.grid(axis="x", alpha=0.3)

    for i, r in enumerate(valid):
        le_str = f"  LE={r.le_adfr:.3f}" if r.le_adfr else ""
        ax.text(r.affinity - 0.08, i,
                f"{r.affinity:.2f}{le_str}", va="center", ha="right", fontsize=7)

    legend_items = [
        mpatches.Patch(color="#2ecc71", label="LE ≥ 0.4 ◎ 優秀"),
        mpatches.Patch(color="#3498db", label="LE ≥ 0.3 ○ 良好"),
        mpatches.Patch(color="#f39c12", label="LE ≥ 0.2 △ 許容"),
        mpatches.Patch(color="#e74c3c", label="LE < 0.2 ✕"),
    ]
    ax.legend(handles=legend_items, loc="lower right", fontsize=8)

    plt.tight_layout()
    out_path = out_dir / "adfr_docking_scores.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  ADFR スコアグラフ保存: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# グラフ 2: 同一分子 ADFR vs smina 並列比較
# ──────────────────────────────────────────────────────────────────────────────

def plot_adfr_vs_smina(results: List[ADFRResult], out_dir: Path):
    """
    同一分子に対して ADFR スコアと smina スコアを左右パネルで比較する。

    これにより「静電項・脱溶媒和の有無でスコアがどう変わるか」が可視化できる。
    """
    if not HAS_MATPLOTLIB:
        return
    _setup_font()

    valid = sorted(
        [r for r in results
         if r.status == "success" and r.affinity < 0 and r.smina_score is not None],
        key=lambda r: r.affinity,
    )
    if not valid:
        return

    n      = len(valid)
    names  = [r.name[:32] + "…" if len(r.name) > 32 else r.name for r in valid]
    adfr_s = [r.affinity       for r in valid]
    smin_s = [r.smina_score    for r in valid]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, max(6, n * 0.5 + 3)))
    fig.suptitle(
        "同一低分子: ADFR vs smina スコア比較\n"
        "(同じ結合部位 / 同じ分子 / 異なる力場)",
        fontsize=13, fontweight="bold",
    )

    # 左: ADFR
    ax1.barh(range(n), adfr_s, color=[_le_color(r.le_adfr) for r in valid], alpha=0.85)
    ax1.set_yticks(range(n))
    ax1.set_yticklabels(names, fontsize=8)
    ax1.set_xlabel("ADFR Affinity (kcal/mol)\n[AutoDock 力場 / 静電+脱溶媒和あり]", fontsize=10)
    ax1.set_title("ADFR (AutoDockFR)", fontsize=11)
    ax1.invert_yaxis()
    ax1.invert_xaxis()
    ax1.grid(axis="x", alpha=0.3)
    for i, r in enumerate(valid):
        ax1.text(r.affinity - 0.06, i, f"{r.affinity:.2f}",
                 va="center", ha="right", fontsize=7)

    # 右: smina
    ax2.barh(range(n), smin_s, color=[_le_color(r.smina_le) for r in valid], alpha=0.85)
    ax2.set_yticks(range(n))
    ax2.set_yticklabels(names, fontsize=8)
    ax2.set_xlabel("smina Affinity (kcal/mol)\n[Vina スコア / 静電・溶媒なし]", fontsize=10)
    ax2.set_title("smina (AutoDock Vina)", fontsize=11)
    ax2.invert_yaxis()
    ax2.invert_xaxis()
    ax2.grid(axis="x", alpha=0.3)
    for i, r in enumerate(valid):
        ax2.text((r.smina_score or 0) - 0.06, i, f"{r.smina_score:.2f}",
                 va="center", ha="right", fontsize=7)

    legend_items = [
        mpatches.Patch(color="#2ecc71", label="LE ≥ 0.4 ◎"),
        mpatches.Patch(color="#3498db", label="LE ≥ 0.3 ○"),
        mpatches.Patch(color="#f39c12", label="LE ≥ 0.2 △"),
        mpatches.Patch(color="#e74c3c", label="LE < 0.2 ✕"),
    ]
    ax1.legend(handles=legend_items, loc="lower right", fontsize=8)
    ax2.legend(handles=legend_items, loc="lower right", fontsize=8)

    fig.text(
        0.5, 0.01,
        "ADFR: 静電相互作用 + 脱溶媒和項を含む AutoDock 力場\n"
        "smina: 経験的 Vina スコア (静電なし / 溶媒なし)\n"
        "スケールが異なるため直接比較不可。同一メソッド内でのランキング比較が有効。",
        ha="center", fontsize=9, color="darkred",
        bbox=dict(boxstyle="round", facecolor="#fff3cd", alpha=0.8),
    )
    plt.tight_layout(rect=[0, 0.07, 1, 1])
    out_path = out_dir / "adfr_vs_smina.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  ADFR vs smina 比較グラフ保存: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# グラフ 3: ADFR (低分子) vs ADCP (環状ペプチド)
# ──────────────────────────────────────────────────────────────────────────────

def plot_adfr_vs_adcp(results: List[ADFRResult], adcp: Optional[dict], out_dir: Path):
    """
    ADFR 低分子スコアと ADCP 環状ペプチドスコアを並列比較。

    両者は同じ AutoDock 力場ベースだが、ADCP は implicit 溶媒を追加するため
    スコアスケールが異なる (ADCP の方が典型的により負)。
    """
    if not HAS_MATPLOTLIB:
        return
    _setup_font()

    valid = sorted(
        [r for r in results if r.status == "success" and r.affinity < 0],
        key=lambda r: r.affinity,
    )
    if not valid:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, max(8, len(valid) * 0.5 + 3)))
    fig.suptitle(
        "AutoDock 力場ベース比較: ADCP (環状ペプチド) vs ADFR (低分子)\n"
        "⚠ 同じ力場ベースだが implicit 溶媒の有無でスコアスケールが異なる",
        fontsize=13, fontweight="bold",
    )

    # ─── 左パネル: ADCP ───
    if adcp and adcp.get("affinity_kcal_mol") and adcp["affinity_kcal_mol"] < 0:
        aff  = adcp["affinity_kcal_mol"]
        seq  = adcp.get("sequence", "GEVDGWATPD")
        le   = adcp.get("le_adcp")
        hac  = adcp.get("hac")

        bar = ax1.barh([f"cyclo-{seq}"], [aff], color="steelblue", alpha=0.85, height=0.4)
        ax1.bar_label(bar, fmt="%.2f kcal/mol", padding=5, fontsize=11, fontweight="bold")

        le_str = f"{le:.4f}" if le else "N/A"
        info = (
            f"配列: cyclo-{seq}\n"
            f"Affinity: {aff:.2f} kcal/mol\n"
            f"HAC: {hac or 'N/A'}  LE: {le_str}\n"
            f"[MC + GB/SA implicit 溶媒]"
        )
        ax1.text(0.03, 0.15, info, transform=ax1.transAxes, fontsize=9, va="bottom",
                 bbox=dict(boxstyle="round", facecolor="lightcyan", alpha=0.6))
        ax1.set_xlim(min(aff * 1.3, -40), 2)
    else:
        ax1.text(0.5, 0.5,
                 "ADCP 結果なし\n→ dock_cyclic_adcp.py を先に実行してください",
                 transform=ax1.transAxes, ha="center", va="center",
                 fontsize=11, color="gray")

    ax1.set_xlabel("ADCP Affinity (kcal/mol)\n[AutoDock 力場 + GB/SA implicit 溶媒]", fontsize=10)
    ax1.set_title("ADCP CrankPep\n(環状ペプチド / implicit 溶媒あり)", fontsize=11)
    ax1.grid(axis="x", alpha=0.3)

    # ─── 右パネル: ADFR ───
    top_n  = min(15, len(valid))
    data   = valid[:top_n]
    names  = [r.name[:40] + "…" if len(r.name) > 40 else r.name for r in data]
    scores = [r.affinity for r in data]
    colors = [_le_color(r.le_adfr) for r in data]

    bars = ax2.barh(range(top_n), scores, color=colors, alpha=0.85)
    ax2.set_yticks(range(top_n))
    ax2.set_yticklabels(names, fontsize=8)
    for i, (r, bar) in enumerate(zip(data, bars)):
        le_str = f"  LE={r.le_adfr:.3f}" if r.le_adfr else ""
        ax2.text(r.affinity - 0.08, i, f"{r.affinity:.2f}{le_str}",
                 va="center", ha="right", fontsize=7)

    ax2.invert_yaxis()
    ax2.invert_xaxis()
    ax2.set_xlabel("ADFR Affinity (kcal/mol)\n[AutoDock 力場 / implicit 溶媒なし]", fontsize=10)
    ax2.set_title(f"ADFR (低分子 Top {top_n})\n(implicit 溶媒なし)", fontsize=11)
    ax2.grid(axis="x", alpha=0.3)

    legend_items = [
        mpatches.Patch(color="#2ecc71", label="LE ≥ 0.4 ◎"),
        mpatches.Patch(color="#3498db", label="LE ≥ 0.3 ○"),
        mpatches.Patch(color="#f39c12", label="LE ≥ 0.2 △"),
        mpatches.Patch(color="#e74c3c", label="LE < 0.2 ✕"),
    ]
    ax2.legend(handles=legend_items, loc="lower right", fontsize=8)

    fig.text(
        0.5, 0.01,
        "ADFR と ADCP は同じ AutoDock 力場を使用 (静電+脱溶媒和+vdW+HBond)。\n"
        "ADCP は GB/SA implicit 溶媒を加えるため典型スコアが低くなる(負に大きい)。\n"
        "LE は同一メソッド内での相対比較に使用してください。",
        ha="center", fontsize=9, color="darkred",
        bbox=dict(boxstyle="round", facecolor="#fff3cd", alpha=0.8),
    )
    plt.tight_layout(rect=[0, 0.08, 1, 1])
    out_path = out_dir / "adfr_vs_adcp_comparison.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  ADFR vs ADCP 比較グラフ保存: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# グラフ 4: LE 3 メソッド比較 (ADCP / ADFR / smina)
# ──────────────────────────────────────────────────────────────────────────────

def plot_le_comparison(results: List[ADFRResult], adcp: Optional[dict], out_dir: Path):
    """
    LE (Ligand Efficiency) を ADCP / ADFR / smina の 3 系統で比較する。

    ADFR と smina の LE は「同じ分子を異なる力場で評価した場合の効率」を示す。
    ADCP の LE は環状ペプチドの効率 (別スケール)。
    """
    if not HAS_MATPLOTLIB:
        return
    _setup_font()

    entries = []

    # ADCP (環状ペプチド)
    if adcp and adcp.get("le_adcp"):
        seq = adcp.get("sequence", "GEVDGWATPD")
        entries.append({
            "label":  f"cyclo-{seq}\n[ADCP]",
            "le":     adcp["le_adcp"],
            "color":  "steelblue",
            "alpha":  0.9,
        })

    # ADFR + smina (同一分子を並べる)
    valid = sorted(
        [r for r in results if r.status == "success" and r.le_adfr is not None],
        key=lambda r: r.affinity,
    )
    for r in valid:
        short = r.name[:26] + "…" if len(r.name) > 26 else r.name
        entries.append({
            "label":  f"{short}\n[ADFR]",
            "le":     r.le_adfr,
            "color":  _le_color(r.le_adfr),
            "alpha":  0.9,
        })
        if r.smina_le is not None:
            entries.append({
                "label":  f"{short}\n[smina]",
                "le":     r.smina_le,
                "color":  _le_color(r.smina_le),
                "alpha":  0.5,   # 半透明で smina を重ねる
            })

    if not entries:
        return

    fig, ax = plt.subplots(figsize=(14, max(6, len(entries) * 0.38 + 2)))
    for i, e in enumerate(entries):
        ax.barh(i, e["le"], color=e["color"], alpha=e["alpha"])
        ax.text(e["le"] + 0.004, i, f"{e['le']:.4f}", va="center", fontsize=7)

    ax.set_yticks(range(len(entries)))
    ax.set_yticklabels([e["label"] for e in entries], fontsize=7)
    ax.set_xlabel("Ligand Efficiency  =  |Affinity| / HAC", fontsize=11)
    ax.set_title(
        "LE 比較: ADCP / ADFR / smina  (3 メソッド)\n"
        "同一メソッド内でのランキング比較が有効",
        fontsize=12,
    )

    ax.axvline(x=0.4, color="green",  linestyle="--", alpha=0.7, label="LE=0.4 ◎ 優秀")
    ax.axvline(x=0.3, color="blue",   linestyle=":",  alpha=0.7, label="LE=0.3 ○ 良好")
    ax.axvline(x=0.2, color="orange", linestyle="-.", alpha=0.7, label="LE=0.2 △ 許容")
    ax.legend(fontsize=9)
    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)

    plt.tight_layout()
    out_path = out_dir / "le_comparison_3methods.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  LE 比較グラフ (3 メソッド) 保存: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# CSV / JSON 出力
# ──────────────────────────────────────────────────────────────────────────────

def write_csv(results: List[ADFRResult], out_dir: Path):
    """ドッキング結果を CSV に保存"""
    csv_path = out_dir / "adfr_docking_results.csv"
    fields = [
        "rank", "name", "adfr_score", "smina_score",
        "hac", "le_adfr", "smina_le",
        "n_clusters", "status", "pair", "dlg_file",
    ]

    # 成功 → スコア順 / 失敗 → 後ろに追加
    success = sorted(
        [r for r in results if r.status == "success" and r.affinity < float("inf")],
        key=lambda r: r.affinity,
    )
    others  = [r for r in results if r not in success]

    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for rank, r in enumerate(success + others, 1):
            w.writerow({
                "rank":        rank,
                "name":        r.name,
                "adfr_score":  r.affinity if r.affinity < float("inf") else "",
                "smina_score": r.smina_score if r.smina_score is not None else "",
                "hac":         r.hac or "",
                "le_adfr":     r.le_adfr or "",
                "smina_le":    r.smina_le or "",
                "n_clusters":  r.n_clusters,
                "status":      r.status,
                "pair":        r.pair,
                "dlg_file":    r.dlg_file or "",
            })
    print(f"  CSV 保存: {csv_path}")


def write_report(results: List[ADFRResult], adcp: Optional[dict], out_dir: Path):
    """JSON レポートを保存"""
    success = sorted(
        [r for r in results if r.status == "success" and r.affinity < float("inf")],
        key=lambda r: r.affinity,
    )
    failed  = [r for r in results if r.status != "success"]

    report = {
        "adfr_small_molecules": [
            {
                "rank":         i + 1,
                "name":         r.name,
                "adfr_score":   r.affinity,
                "smina_score":  r.smina_score,
                "hac":          r.hac,
                "le_adfr":      r.le_adfr,
                "smina_le":     r.smina_le,
                "n_clusters":   r.n_clusters,
                "pair":         r.pair,
                "dlg_file":     r.dlg_file,
                "pose_file":    r.pose_file,
            }
            for i, r in enumerate(success)
        ],
        "adcp_cyclic_peptide": adcp,
        "failed_molecules": [
            {"name": r.name, "status": r.status, "error": r.error_message}
            for r in failed
        ],
        "scoring_info": {
            "ADFR": {
                "tool":          "AutoDockFR (adfr v0.0.22)",
                "method":        "Genetic Algorithm",
                "force_field":   "AutoDock 4 (静電 + 脱溶媒和 + vdW + HBond)",
                "solvent":       "なし (真空中)",
                "typical_range": "-5 〜 -15 kcal/mol",
            },
            "ADCP": {
                "tool":          "AutoDock CrankPep (adcp v0.0.25)",
                "method":        "Monte Carlo (-cyc, 環状ペプチド専用)",
                "force_field":   "AutoDock 4 + OpenMM GB/SA implicit 溶媒",
                "solvent":       "GB/SA implicit 溶媒",
                "typical_range": "-10 〜 -35 kcal/mol",
            },
            "smina": {
                "tool":          "smina (AutoDock Vina based)",
                "method":        "Vina gradient descent",
                "force_field":   "経験的 Vina スコア (静電なし / 溶媒なし)",
                "solvent":       "なし",
                "typical_range": "-5 〜 -10 kcal/mol",
            },
            "comparison_note": (
                "ADFR と ADCP は同じ AutoDock 力場ベース (静電+脱溶媒和+vdW+HBond)。"
                "ADCP は GB/SA implicit 溶媒を加えるためスコアが低くなる傾向がある。"
                "smina(Vina) は静電・溶媒なしの経験的スコア。"
                "異なるメソッド間の直接スコア比較は不可。"
                "LE (Ligand Efficiency) も同一メソッド内での相対比較に使用してください。"
            ),
        },
    }

    out_path = out_dir / "adfr_comparison_report.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"  レポート保存: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# コンソールサマリー
# ──────────────────────────────────────────────────────────────────────────────

def print_summary(results: List[ADFRResult], adcp: Optional[dict]):
    success = sorted(
        [r for r in results if r.status == "success" and r.affinity < float("inf")],
        key=lambda r: r.affinity,
    )
    failed  = [r for r in results if r.status != "success"]

    print()
    print("=" * 96)
    print("  ADFR ドッキング結果  [AutoDock 力場 / 静電+脱溶媒和あり / implicit 溶媒なし]")
    print("=" * 96)
    print(f"  {'順位':<4}  {'ADFR':>10}  {'smina':>8}  {'HAC':>5}  "
          f"{'LE(ADFR)':>9}  {'LE(smina)':>9}  分子名")
    print("  " + "─" * 94)
    for i, r in enumerate(success, 1):
        smina_str = f"{r.smina_score:>8.3f}" if r.smina_score is not None else "     N/A"
        le_a_str  = f"{r.le_adfr:.4f}"       if r.le_adfr     is not None else "    N/A"
        le_s_str  = f"{r.smina_le:.4f}"      if r.smina_le    is not None else "    N/A"
        print(f"  {i:<4}  {r.affinity:>10.3f}  {smina_str}  "
              f"{str(r.hac or '?'):>5}  {le_a_str:>9}  {le_s_str:>9}  {r.name[:42]}")

    if failed:
        print(f"\n  失敗 ({len(failed)} 件):")
        for r in failed:
            print(f"    ✗ {r.name[:50]}: {r.status}  {r.error_message or ''}")

    if adcp and adcp.get("affinity_kcal_mol"):
        seq  = adcp.get("sequence", "GEVDGWATPD")
        aff  = adcp["affinity_kcal_mol"]
        le   = adcp.get("le_adcp")
        print()
        print("  ─" * 48)
        print(f"  ADCP 環状ペプチド: cyclo-{seq}")
        print(f"    Affinity: {aff:.3f} kcal/mol  [implicit 溶媒あり]")
        print(f"    LE(ADCP): {le:.4f}" if le else "    LE(ADCP): N/A")

    print()
    print("  ─" * 48)
    print("  スコアスケールの目安:")
    print("    ADFR : AutoDock 力場 (静電+脱溶媒和)  / implicit 溶媒なし  → -5〜-15  kcal/mol")
    print("    ADCP : AutoDock 力場 + GB/SA implicit 溶媒                  → -10〜-35 kcal/mol")
    print("    smina: Vina スコア  (静電・溶媒なし)                          → -5〜-10  kcal/mol")
    print("  異なるメソッド間のスコア直接比較は不可。LE も同一メソッド内での使用を推奨。")
    print("  ─" * 48)


# ──────────────────────────────────────────────────────────────────────────────
# メイン
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "AutoDockFR (ADFR) で低分子をドッキングし、"
            "ADCP 環状ペプチド結果と比較する。"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--top-n", type=int, default=15,
        help="ドッキングする分子数 (デフォルト: 15)",
    )
    parser.add_argument(
        "--quick", action="store_true",
        help=f"クイックモード: N={QUICK_N_RUNS}, n={QUICK_N_EVALS} (精度は低い)",
    )
    parser.add_argument(
        "--n-runs", type=int, default=DEFAULT_N_RUNS,
        help=f"GA 探索回数 (デフォルト: {DEFAULT_N_RUNS})",
    )
    parser.add_argument(
        "--n-evals", type=int, default=DEFAULT_N_EVALS,
        help=f"評価ステップ数 (デフォルト: {DEFAULT_N_EVALS})",
    )
    parser.add_argument(
        "--timeout", type=int, default=600,
        help="1 分子あたりのタイムアウト秒数 (デフォルト: 600)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="既存結果を上書きして再実行する",
    )
    args = parser.parse_args()

    if args.quick:
        args.n_runs  = QUICK_N_RUNS
        args.n_evals = QUICK_N_EVALS
        args.timeout = 180
        print(f"  [クイックモード] N={args.n_runs}, n={args.n_evals}, timeout={args.timeout}s")

    print("\n" + "=" * 65)
    print("  AutoDockFR (ADFR) 低分子ドッキング + ADCP 比較")
    print("=" * 65)

    # ① micromamba 確認
    mamba = find_micromamba()
    if mamba is None:
        print("[エラー] micromamba が見つかりません。")
        print("  インストール: curl -Ls https://micro.mamba.pm/install.sh | bash")
        sys.exit(1)
    print(f"  micromamba : {mamba}")

    # ② adcpsuite 環境確認
    if not check_adcpsuite_env(mamba):
        print("[エラー] adcpsuite 環境が見つかりません。")
        print("  セットアップ: bash adcpsuite_micromamba.sh")
        sys.exit(1)
    print("  adcpsuite  : OK")

    ADFR_DIR.mkdir(parents=True, exist_ok=True)
    pdbqt_dir = ADFR_DIR / "ligands_pdbqt"
    pdbqt_dir.mkdir(exist_ok=True)

    # ③ receptor.trg 確認
    trg_file = find_receptor_trg()
    if trg_file is None:
        print("[エラー] receptor.trg が見つかりません。")
        print("  実行: python dock_cyclic_adcp.py --setup-receptor")
        sys.exit(1)
    size_mb = trg_file.stat().st_size / 1024 / 1024
    print(f"  receptor.trg : {trg_file}  ({size_mb:.1f} MB)")

    # ④ 分子読み込み
    print(f"\n  ① 低分子候補読み込み中 (Top {args.top_n}) ...")
    try:
        molecules = load_molecules(top_n=args.top_n)
    except FileNotFoundError as e:
        print(f"[エラー] {e}")
        sys.exit(1)
    print(f"     {len(molecules)} 件ロード完了")

    # ⑤ ADFR ドッキング (分子ごとループ)
    print(f"\n  ② ADFR ドッキング開始")
    print(f"     分子数: {len(molecules)} 件  /  N={args.n_runs}  /  n={args.n_evals}  /  timeout={args.timeout}s")

    results: List[ADFRResult] = []
    for i, mol in enumerate(molecules, 1):
        name     = mol["name"]
        sdf_path = mol["sdf_path"]
        print(f"\n  [{i:02d}/{len(molecules):02d}] {name[:58]}")

        # SDF → PDBQT
        safe  = re.sub(r"[^\w\-]", "_", name)[:50].lower()
        pdbqt = pdbqt_dir / f"{safe}.pdbqt"
        print(f"    変換: SDF → PDBQT ...")
        ok = convert_sdf_to_pdbqt(mamba, sdf_path, pdbqt, force=args.force)
        if not ok:
            print(f"    ✗ PDBQT 変換失敗 ({sdf_path.name})")
            results.append(ADFRResult(
                name=name, affinity=float("inf"), best_energy=float("inf"),
                n_clusters=0, status="convert_failed",
                error_message="SDF→PDBQT 変換失敗",
                hac=mol.get("hac"), smina_score=mol.get("smina_score"),
                smina_le=mol.get("smina_le"), pair=mol.get("pair", ""),
            ))
            continue

        # ADFR ドッキング
        print(f"    ADFR ドッキング中 ...")
        result = run_adfr_docking(
            mamba=mamba, name=name, pdbqt_path=pdbqt, trg_file=trg_file,
            n_runs=args.n_runs, n_evals=args.n_evals,
            timeout=args.timeout, force=args.force,
            hac=mol.get("hac"), smina_score=mol.get("smina_score"),
            smina_le=mol.get("smina_le"), pair=mol.get("pair", ""),
        )
        results.append(result)

        if result.status == "success" and result.affinity < float("inf"):
            le_str = f"  LE={result.le_adfr:.4f}" if result.le_adfr else ""
            print(f"    ✓ Affinity: {result.affinity:.3f} kcal/mol{le_str}")
        else:
            print(f"    ✗ {result.status}: {result.error_message or ''}")

    # ⑥ ADCP 結果読み込み
    print("\n  ③ ADCP 結果読み込み中 ...")
    adcp = load_adcp_result()
    if adcp and adcp.get("affinity_kcal_mol"):
        seq = adcp.get("sequence", "?")
        aff = adcp["affinity_kcal_mol"]
        print(f"     ADCP: cyclo-{seq}  →  {aff:.3f} kcal/mol  [implicit 溶媒あり]")
    else:
        print("     ADCP 結果なし (dock_cyclic_adcp.py を先に実行することを推奨)")

    # ⑦ グラフ生成
    print("\n  ④ グラフ生成中 ...")
    plot_adfr_scores(results, ADFR_DIR)
    plot_adfr_vs_smina(results, ADFR_DIR)
    plot_adfr_vs_adcp(results, adcp, ADFR_DIR)
    plot_le_comparison(results, adcp, ADFR_DIR)

    # ⑧ CSV / JSON 保存
    write_csv(results, ADFR_DIR)
    write_report(results, adcp, ADFR_DIR)

    # ⑨ サマリー表示
    print_summary(results, adcp)

    n_ok   = sum(1 for r in results if r.status == "success" and r.affinity < float("inf"))
    n_fail = len(results) - n_ok
    print(f"\n  完了!  成功: {n_ok} 件  /  失敗: {n_fail} 件  /  出力先: {ADFR_DIR}")


if __name__ == "__main__":
    main()
