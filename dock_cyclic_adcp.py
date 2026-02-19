#!/usr/bin/env python3
"""
dock_cyclic_adcp.py
===================
AutoDock CrankPep (ADCP) を使って環状ペプチド cyclo-GEVDGWATPD をドッキングし、
smina による低分子ドッキング結果と比較する。

【スコアスケールの違いについて】
ADCP と smina は全く異なるスコアリング関数を使用するため、スコアの直接比較不可:
  - ADCP  : Monte Carlo + implicit 溶媒 (環状ペプチド専用) → 典型: -10〜-35 kcal/mol
  - smina : AutoDock Vina スコア (低分子)                   → 典型: -5〜-10 kcal/mol
LE (Ligand Efficiency = |Affinity| / HAC) も参考値として表示しますが、
scoring function が異なるため厳密な比較には使えません。

【必要環境】
  micromamba + adcpsuite 環境 (ADCP v0.0.25, AGFR v0.0.22)
  セットアップ: bash adcpsuite_micromamba.sh

Usage:
  # 環境セットアップ (初回のみ)
  bash adcpsuite_micromamba.sh

  # ターゲットファイル作成 (初回のみ)
  python dock_cyclic_adcp.py --setup-receptor

  # 通常ドッキング (N=5, n=100000, ≈5〜10分)
  python dock_cyclic_adcp.py

  # クイックテスト (N=1, n=5000, ≈30秒〜1分)
  python dock_cyclic_adcp.py --quick

  # 既存結果を強制上書き
  python dock_cyclic_adcp.py --force

ADCPコマンド参考:
  adcp -O -T receptor.trg -s GEVDGWATPD -cyc -N 5 -n 100000
       -nmin 5 -nitr 500 -env implicit -dr -reint
       -o result_gevdgwatpd -w results/adcp_docking/gevdgwatpd/
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
from dataclasses import dataclass
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
DOCKING_DIR = RESULTS_DIR / "docking"
BRIDGE_DIR  = RESULTS_DIR / "bridge"
ADCP_DIR    = RESULTS_DIR / "adcp_docking"

DEFAULT_SEQUENCE = "GEVDGWATPD"

# cyclo-GEVDGWATPD の重原子数 (compare_cyclic_peptide.py で計算済み)
CYCLIC_GEVDGWATPD_HAC = 73

# PRODIGY 結合自由エネルギー (参考値)
PRODIGY_DG = -14.7  # kcal/mol


# ──────────────────────────────────────────────────────────────────────────────
# データクラス
# ──────────────────────────────────────────────────────────────────────────────

@dataclass
class ADCPResult:
    """ADCP ドッキング結果"""
    sequence:      str
    affinity:      float          # Best affinity (kcal/mol) from mode 1
    best_energy:   float          # Best MC energy
    n_clusters:    int
    status:        str            # success / failed / timeout / error
    error_message: Optional[str] = None
    dlg_file:      Optional[str] = None
    pdb_file:      Optional[str] = None
    hac:           Optional[int] = None
    le:            Optional[float] = None  # |affinity| / HAC


# ──────────────────────────────────────────────────────────────────────────────
# 環境チェック
# ──────────────────────────────────────────────────────────────────────────────

def find_micromamba() -> Optional[str]:
    """micromamba の実行パスを探す"""
    # 1. which コマンド
    try:
        r = subprocess.run(["which", "micromamba"], capture_output=True, text=True)
        if r.returncode == 0:
            return r.stdout.strip()
    except Exception:
        pass

    # 2. よくあるパス
    candidates = [
        os.path.expanduser("~/.local/bin/micromamba"),
        os.path.expanduser("~/micromamba/bin/micromamba"),
        os.path.expanduser("~/mambaforge/bin/micromamba"),
        "/usr/local/bin/micromamba",
        "/opt/homebrew/bin/micromamba",
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    return None


def find_adcpsuite_env(mamba: str) -> Optional[str]:
    """adcpsuite 環境のフルパスを返す（異なるconda/mambaルートにも対応）"""
    # 1. micromamba env list のパス列から探す
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

    # 2. よくあるインストール先を直接探す
    candidates = [
        os.path.expanduser("~/miniforge3/envs/adcpsuite"),
        os.path.expanduser("~/miniforge/envs/adcpsuite"),
        os.path.expanduser("~/mambaforge/envs/adcpsuite"),
        os.path.expanduser("~/micromamba/envs/adcpsuite"),
        os.path.expanduser("~/miniconda3/envs/adcpsuite"),
        os.path.expanduser("~/opt/anaconda3/envs/adcpsuite"),
        "/opt/anaconda3/envs/adcpsuite",
        "/opt/homebrew/Caskroom/miniforge/base/envs/adcpsuite",
    ]
    for p in candidates:
        if os.path.isdir(p):
            return p
    return None


def check_adcpsuite_env(mamba: str) -> bool:
    """adcpsuite conda 環境が存在するか確認"""
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
        full_cmd,
        shell=True,
        capture_output=True,
        text=True,
        cwd=cwd,
        timeout=timeout,
    )


# ──────────────────────────────────────────────────────────────────────────────
# 受容体準備: receptor.pdb → receptor.trg
# ──────────────────────────────────────────────────────────────────────────────

def setup_receptor(mamba: str, force: bool = False) -> Path:
    """
    receptor.pdb から ADCP ターゲットファイル (.trg) を作成する。

    手順:
      1. receptor.pdb  → receptor_rec.pdbqt    (agfr --toPdbqt)
      2. peptide_ref.pdb → peptide_ref.pdbqt   (obabel)
      3. receptor_rec.pdbqt + peptide_ref.pdbqt → receptor.trg  (agfr)

    Parameters
    ----------
    force : bool
        True の場合、既存の .trg ファイルを再作成する。
    """
    ADCP_DIR.mkdir(parents=True, exist_ok=True)
    trg_file = ADCP_DIR / "receptor.trg"

    if trg_file.exists() and not force:
        size_mb = trg_file.stat().st_size / 1024 / 1024
        print(f"  ターゲットファイル既存 ({size_mb:.1f} MB): {trg_file}")
        return trg_file

    # 入力ファイル確認
    receptor_pdb = DOCKING_DIR / "receptor.pdb"
    peptide_pdb  = DOCKING_DIR / "peptide_ref.pdb"

    if not receptor_pdb.exists():
        raise FileNotFoundError(
            f"receptor.pdb が見つかりません: {receptor_pdb}\n"
            f"先に smina ドッキング (dock_with_smina.py) を実行してください。"
        )
    if not peptide_pdb.exists():
        raise FileNotFoundError(
            f"peptide_ref.pdb が見つかりません: {peptide_pdb}\n"
            f"先に smina ドッキング (dock_with_smina.py) を実行してください。"
        )

    # 作業ディレクトリにコピー (相対パス問題を避けるため)
    shutil.copy2(receptor_pdb, ADCP_DIR / "receptor.pdb")
    shutil.copy2(peptide_pdb,  ADCP_DIR / "peptide_ref.pdb")

    cwd = str(ADCP_DIR)

    # ── Step 1: receptor.pdb → receptor_rec.pdbqt ──
    print("  Step 1/3: receptor.pdb → receptor_rec.pdbqt  (agfr --toPdbqt)")
    r = run_in_adcpsuite(mamba, "agfr -r receptor.pdb --toPdbqt", cwd=cwd, timeout=300)
    if r.returncode != 0:
        print(f"  STDOUT: {r.stdout[:500]}")
        raise RuntimeError(f"agfr --toPdbqt 失敗:\n{r.stderr[:1000]}")

    rec_pdbqt = ADCP_DIR / "receptor_rec.pdbqt"
    if not rec_pdbqt.exists():
        raise FileNotFoundError("receptor_rec.pdbqt が生成されませんでした。agfr のログを確認してください。")

    # ── Step 2: peptide_ref.pdb → peptide_ref.pdbqt ──
    print("  Step 2/3: peptide_ref.pdb → peptide_ref.pdbqt  (obabel)")
    r = run_in_adcpsuite(
        mamba,
        "obabel -ipdb peptide_ref.pdb -opdbqt -O peptide_ref.pdbqt --partialcharge gasteiger",
        cwd=cwd, timeout=120,
    )
    if r.returncode != 0 or not (ADCP_DIR / "peptide_ref.pdbqt").exists():
        # obabel 失敗時は prepare_ligand4.py を試みる
        print("    obabel 失敗。prepare_ligand4.py を試みます...")
        r = run_in_adcpsuite(
            mamba,
            "prepare_ligand4.py -l peptide_ref.pdb -o peptide_ref.pdbqt",
            cwd=cwd, timeout=120,
        )
        if r.returncode != 0 or not (ADCP_DIR / "peptide_ref.pdbqt").exists():
            raise RuntimeError(f"peptide_ref.pdb → PDBQT 変換失敗:\n{r.stderr[:1000]}")

    # ── Step 3: receptor_rec.pdbqt + peptide_ref.pdbqt → receptor.trg ──
    print("  Step 3/3: receptor_rec.pdbqt + peptide_ref.pdbqt → receptor.trg  (agfr)")
    r = run_in_adcpsuite(
        mamba,
        "agfr -r receptor_rec.pdbqt -l peptide_ref.pdbqt -asv 1.1 -o receptor",
        cwd=cwd, timeout=600,
    )
    if r.returncode != 0:
        print(f"  STDOUT: {r.stdout[:500]}")
        raise RuntimeError(f"agfr ターゲット作成失敗:\n{r.stderr[:1000]}")

    if not trg_file.exists():
        raise FileNotFoundError("receptor.trg が生成されませんでした。agfr のログを確認してください。")

    size_mb = trg_file.stat().st_size / 1024 / 1024
    print(f"  ターゲットファイル作成完了 ({size_mb:.1f} MB): {trg_file}")
    return trg_file


# ──────────────────────────────────────────────────────────────────────────────
# ADCP ドッキング
# ──────────────────────────────────────────────────────────────────────────────

def run_adcp_docking(
    mamba:    str,
    sequence: str,
    trg_file: Path,
    n_runs:   int = 5,
    n_evals:  int = 100000,
    n_min:    int = 5,
    n_itr:    int = 500,
    timeout:  int = 3600,
    force:    bool = False,
) -> ADCPResult:
    """
    ADCP で環状ペプチドをドッキングする。

    Parameters
    ----------
    sequence : str
        アミノ酸配列 (1文字コード, 例: "GEVDGWATPD")。大文字・小文字どちらでも可。
    trg_file : Path
        受容体ターゲットファイル (.trg)。setup_receptor() で作成。
    n_runs   : int  MC 探索の実行回数 (デフォルト: 5)
    n_evals  : int  各探索の評価ステップ数 (デフォルト: 100000)
    n_min    : int  OpenMM 最小化回数 (デフォルト: 5)
    n_itr    : int  ADCP イテレーション数 (デフォルト: 500)
    timeout  : int  タイムアウト秒数 (デフォルト: 3600)
    force    : bool 既存結果を上書きするか
    """
    seq_upper = sequence.upper()
    seq_lower = sequence.lower()
    job_name  = seq_lower

    out_dir  = ADCP_DIR / job_name
    out_dir.mkdir(parents=True, exist_ok=True)

    dlg_file = out_dir / f"result_{job_name}_summary.dlg"
    pdb_file = out_dir / f"result_{job_name}_out.pdb"

    # キャッシュ確認
    if dlg_file.exists() and not force:
        print(f"  既存ドッキング結果を使用: {dlg_file}")
        affinity, best_energy, n_clusters = _parse_dlg_file(dlg_file)
        hac = _calc_cyclic_hac(seq_upper)
        le  = round(abs(affinity) / hac, 4) if (hac and affinity < 0) else None
        return ADCPResult(
            sequence=seq_upper, affinity=affinity, best_energy=best_energy,
            n_clusters=n_clusters, status="success",
            dlg_file=str(dlg_file), pdb_file=str(pdb_file) if pdb_file.exists() else None,
            hac=hac, le=le,
        )

    # ADCP コマンド構築
    # 参考: adcp -O -T receptor.trg -s GEVDGWATPD -cyc
    #            -N 5 -n 100000 -nmin 5 -nitr 500 -env implicit -dr -reint
    #            -o result_gevdgwatpd -w <output_dir>
    cmd = (
        f'adcp -O'
        f' -T "{trg_file.resolve()}"'
        f' -s "{seq_upper}"'
        f' -cyc'
        f' -N {n_runs}'
        f' -n {n_evals}'
        f' -nmin {n_min}'
        f' -nitr {n_itr}'
        f' -env implicit'
        f' -dr'
        f' -reint'
        f' -o result_{job_name}'
        f' -w "{out_dir.resolve()}"'
    )

    print(f"  コマンド: micromamba run -n adcpsuite {cmd}")

    try:
        r = run_in_adcpsuite(mamba, cmd, timeout=timeout)

        if r.returncode != 0:
            return ADCPResult(
                sequence=seq_upper, affinity=float("inf"),
                best_energy=float("inf"), n_clusters=0,
                status="failed", error_message=r.stderr[:500],
            )

        affinity, best_energy, n_clusters = _parse_dlg_file(dlg_file)
        hac = _calc_cyclic_hac(seq_upper)
        le  = round(abs(affinity) / hac, 4) if (hac and affinity < 0) else None

        return ADCPResult(
            sequence=seq_upper, affinity=affinity, best_energy=best_energy,
            n_clusters=n_clusters, status="success",
            dlg_file=str(dlg_file), pdb_file=str(pdb_file) if pdb_file.exists() else None,
            hac=hac, le=le,
        )

    except subprocess.TimeoutExpired:
        return ADCPResult(
            sequence=seq_upper, affinity=float("inf"),
            best_energy=float("inf"), n_clusters=0,
            status="timeout", error_message=f"Timeout after {timeout}s",
        )
    except Exception as e:
        return ADCPResult(
            sequence=seq_upper, affinity=float("inf"),
            best_energy=float("inf"), n_clusters=0,
            status="error", error_message=str(e),
        )


def _parse_dlg_file(dlg_file: Path) -> Tuple[float, float, int]:
    """
    ADCP の .dlg ファイルから affinity, best_energy, n_clusters を抽出する。

    .dlg フォーマット例:
        mode |  affinity  | ...
             | (kcal/mol) | ...
        -----+------------+---
           1       -25.3  ...
        bestEnergy in run 2 -26.1 (0)
    """
    affinity    = float("inf")
    best_energy = float("inf")
    n_clusters  = 0

    if not dlg_file.exists():
        return affinity, best_energy, n_clusters

    content = dlg_file.read_text(errors="replace")

    # Best affinity: mode 1 行の最初の数値
    m = re.search(r"^\s+1\s+([-\d\.]+)", content, re.MULTILINE)
    if m:
        affinity = float(m.group(1))

    # Best energy from MC search
    m = re.search(r"bestEnergy in run \d+ ([-\d\.]+)", content)
    if m:
        best_energy = float(m.group(1))

    # Cluster count (テーブルの行数をカウント)
    n_clusters = len(re.findall(r"^\s+\d+\s+[-\d\.]+", content, re.MULTILINE))

    return affinity, best_energy, n_clusters


def _calc_cyclic_hac(sequence: str) -> Optional[int]:
    """環状ペプチドの重原子数 (Heavy Atom Count) を返す"""
    if sequence.upper() == "GEVDGWATPD":
        return CYCLIC_GEVDGWATPD_HAC

    # compare_cyclic_peptide.py の関数を利用
    try:
        from rdkit import Chem
        sys.path.insert(0, str(BASE_DIR))
        from compare_cyclic_peptide import build_cyclic_peptide_smiles
        smiles = build_cyclic_peptide_smiles(sequence)
        mol = Chem.MolFromSmiles(smiles)
        return mol.GetNumAtoms() if mol else None
    except Exception:
        return None


# ──────────────────────────────────────────────────────────────────────────────
# smina 結果読み込み
# ──────────────────────────────────────────────────────────────────────────────

def load_smina_results(top_n: int = 15) -> List[dict]:
    """
    既存の smina ドッキング結果から上位 top_n 件を読み込む。
    Result_Best/summary.csv を優先して使用する。
    """
    # ── Result_Best/summary.csv を優先 ──
    summary_csv = BASE_DIR / "Result_Best" / "summary.csv"
    if summary_csv.exists():
        records = []
        with open(summary_csv) as f:
            for row in csv.DictReader(f):
                try:
                    score = float(row["score"])
                    hac   = int(row["hac"])   if row.get("hac")  else None
                    le    = float(row["le"])   if row.get("le")   else None
                    records.append({
                        "name":   row["name"],
                        "score":  score,
                        "hac":    hac,
                        "le":     le,
                        "pair":   row.get("pair", ""),
                        "method": "smina",
                    })
                except (ValueError, KeyError):
                    pass
        if records:
            records.sort(key=lambda r: r["score"])
            return records[:top_n]

    # ── フォールバック: 個別 CSV を集約 ──
    from rdkit import Chem

    def _hac_from_smiles(smiles: str) -> Optional[int]:
        try:
            mol = Chem.MolFromSmiles(smiles)
            return mol.GetNumAtoms() if mol else None
        except Exception:
            return None

    records = []

    # peptidomimetics (STEP 3)
    docking_csv = DOCKING_DIR / "docking_results.csv"
    if docking_csv.exists():
        with open(docking_csv) as f:
            for row in csv.DictReader(f):
                score_str = row.get("best_score_kcal/mol", "N/A")
                if score_str == "N/A":
                    continue
                try:
                    score = float(score_str)
                    hac   = _hac_from_smiles(row.get("smiles", ""))
                    le    = round(abs(score) / hac, 4) if hac else None
                    records.append({
                        "name": row["name"], "score": score,
                        "hac": hac, "le": le,
                        "pair": "TRP6-ALA7-ASP4", "method": "smina",
                    })
                except ValueError:
                    pass

    # bridge (STEP 7)
    if BRIDGE_DIR.is_dir():
        for csv_path in sorted(BRIDGE_DIR.glob("bridge_*_results.csv")):
            with open(csv_path) as f:
                for row in csv.DictReader(f):
                    score_str = row.get("dock_score", "").strip()
                    if not score_str:
                        continue
                    if row.get("DrugLike", "True").strip() != "True":
                        continue
                    try:
                        score = float(score_str)
                        hac   = _hac_from_smiles(row.get("smiles", ""))
                        le    = round(abs(score) / hac, 4) if hac else None
                        records.append({
                            "name":  row["name"],
                            "score": score, "hac": hac, "le": le,
                            "pair":  f"{row.get('res1','')}-{row.get('res2','')}",
                            "method": "smina",
                        })
                    except ValueError:
                        pass

    # フィルタ・重複除去・ソート
    records = [r for r in records if r["score"] <= -5.0]
    records.sort(key=lambda r: r["score"])
    seen, unique = set(), []
    for r in records:
        if r["name"] not in seen:
            seen.add(r["name"])
            unique.append(r)
    return unique[:top_n]


# ──────────────────────────────────────────────────────────────────────────────
# 可視化
# ──────────────────────────────────────────────────────────────────────────────

def plot_comparison(adcp: ADCPResult, smina: List[dict], out_dir: Path):
    """
    ADCP (環状ペプチド) vs smina (低分子) の比較グラフを 2 パネルで生成する。

    注意: ADCP と smina は異なるスコアリング関数のため y 軸スケールが異なる。
          グラフは参考情報として提供。
    """
    if not HAS_MATPLOTLIB:
        print("  [警告] matplotlib が利用できません")
        return

    try:
        from plot_utils import setup_japanese_font
        setup_japanese_font()
    except ImportError:
        pass

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
    fig.suptitle(
        "環状ペプチド (ADCP CrankPep) vs 低分子 (smina) 結合親和性比較\n"
        "⚠ 異なるスコアリング関数のため数値の直接比較は不可",
        fontsize=13, fontweight="bold",
    )

    # ─── 左パネル: ADCP スコア (環状ペプチド) ───
    if adcp.status == "success" and adcp.affinity < 0:
        bar = ax1.barh(
            [f"cyclo-{adcp.sequence}"],
            [adcp.affinity],
            color="steelblue", alpha=0.85, height=0.4,
        )
        ax1.bar_label(bar, fmt="%.2f kcal/mol", padding=5, fontsize=11, fontweight="bold")

        # PRODIGY ΔG 参考線
        ax1.axvline(
            x=PRODIGY_DG, color="red", linestyle="--", alpha=0.7,
            label=f"PRODIGY ΔG ({PRODIGY_DG} kcal/mol)",
        )
        ax1.legend(fontsize=9, loc="lower right")

        xlim = min(adcp.affinity * 1.3, -40)
        ax1.set_xlim(xlim, 2)

        # 情報テキスト
        le_str = f"{adcp.le:.4f}" if adcp.le else "N/A"
        info = (
            f"配列: cyclo-{adcp.sequence}\n"
            f"Affinity (ADCP): {adcp.affinity:.2f} kcal/mol\n"
            f"HAC: {adcp.hac or 'N/A'}\n"
            f"LE (ADCP): {le_str}\n"
            f"Clusters: {adcp.n_clusters}"
        )
        ax1.text(
            0.03, 0.15, info, transform=ax1.transAxes,
            fontsize=9, va="bottom",
            bbox=dict(boxstyle="round", facecolor="lightcyan", alpha=0.6),
        )
    else:
        ax1.text(
            0.5, 0.5, f"ドッキング失敗\n({adcp.status})\n{adcp.error_message or ''}",
            transform=ax1.transAxes, ha="center", va="center",
            fontsize=11, color="red",
        )

    ax1.set_xlabel("ADCP Affinity (kcal/mol)", fontsize=11)
    ax1.set_title(
        "ADCP CrankPep\n[環状ペプチド専用 / MC + implicit 溶媒]",
        fontsize=11,
    )
    ax1.grid(axis="x", alpha=0.3)

    # ─── 右パネル: smina スコア (低分子 Top 15) ───
    if smina:
        top_n  = min(15, len(smina))
        data   = smina[:top_n]
        names  = [r["name"][:38] + "…" if len(r["name"]) > 38 else r["name"] for r in data]
        scores = [r["score"] for r in data]

        # LE で色分け
        def _le_color(le):
            if le is None:    return "gray"
            if le >= 0.4:     return "#2ecc71"
            if le >= 0.3:     return "#3498db"
            if le >= 0.2:     return "#f39c12"
            return "#e74c3c"

        colors = [_le_color(r.get("le")) for r in data]
        bars   = ax2.barh(range(top_n), scores, color=colors, alpha=0.85)
        ax2.set_yticks(range(top_n))
        ax2.set_yticklabels(names, fontsize=8)

        # ラベル
        for i, (r, bar) in enumerate(zip(data, bars)):
            le_str = f" LE={r['le']:.3f}" if r.get("le") else ""
            ax2.text(
                r["score"] - 0.05, i,
                f"{r['score']:.2f}{le_str}",
                va="center", ha="right", fontsize=7,
            )

        ax2.invert_yaxis()
        ax2.invert_xaxis()

        # LE 凡例
        legend_items = [
            mpatches.Patch(color="#2ecc71", label="LE ≥ 0.4 ◎ 優秀"),
            mpatches.Patch(color="#3498db", label="LE ≥ 0.3 ○ 良好"),
            mpatches.Patch(color="#f39c12", label="LE ≥ 0.2 △ 許容"),
            mpatches.Patch(color="#e74c3c", label="LE < 0.2 ✕"),
        ]
        ax2.legend(handles=legend_items, loc="lower right", fontsize=8)

    ax2.set_xlabel("smina Affinity (kcal/mol)", fontsize=11)
    ax2.set_title(
        f"smina (AutoDock Vina)\n[低分子 Top {min(15, len(smina))} / DrugLike フィルタ済み]",
        fontsize=11,
    )
    ax2.grid(axis="x", alpha=0.3)

    # 注意書き
    fig.text(
        0.5, 0.01,
        "⚠ ADCP (環状ペプチド専用、implicit溶媒) と smina (Vina スコア) は異なるスコアリング関数のため直接比較不可。\n"
        "LE も異なる scoring function 間では厳密な比較には使えません (参考値)。",
        ha="center", fontsize=9, color="darkred",
        bbox=dict(boxstyle="round", facecolor="#fff3cd", alpha=0.8),
    )

    plt.tight_layout(rect=[0, 0.05, 1, 1])
    out_path = out_dir / "adcp_vs_smina_comparison.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  比較グラフ保存: {out_path}")


def plot_le_comparison(adcp: ADCPResult, smina: List[dict], out_dir: Path):
    """LE (Ligand Efficiency) 比較グラフ"""
    if not HAS_MATPLOTLIB:
        return

    try:
        from plot_utils import setup_japanese_font
        setup_japanese_font()
    except ImportError:
        pass

    entries = []

    # ADCP エントリ
    if adcp.status == "success" and adcp.le is not None:
        entries.append({
            "label":  f"cyclo-{adcp.sequence} [ADCP]",
            "le":     adcp.le,
            "color":  "steelblue",
        })

    # smina エントリ
    for r in smina:
        if r.get("le") is not None:
            short = r["name"][:28] + "…" if len(r["name"]) > 28 else r["name"]
            le = r["le"]
            entries.append({
                "label": f"{short} [smina]",
                "le":    le,
                "color": "#2ecc71" if le >= 0.4 else
                         "#3498db" if le >= 0.3 else
                         "#f39c12" if le >= 0.2 else "#e74c3c",
            })

    if not entries:
        print("  [情報] LE グラフ: LE が計算できるデータがありません")
        return

    fig, ax = plt.subplots(figsize=(14, max(5, len(entries) * 0.45 + 2)))

    labels = [e["label"] for e in entries]
    les    = [e["le"]    for e in entries]
    colors = [e["color"] for e in entries]

    bars = ax.barh(range(len(entries)), les, color=colors, alpha=0.85)
    ax.set_yticks(range(len(entries)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel("Ligand Efficiency  =  |Affinity| / HAC", fontsize=11)
    ax.set_title(
        "Ligand Efficiency 比較\n"
        "(ADCP と smina は異なる scoring function のため参考値)",
        fontsize=12,
    )

    # 参考線
    ax.axvline(x=0.4, color="green",  linestyle="--", alpha=0.7, label="LE=0.4 ◎ 優秀")
    ax.axvline(x=0.3, color="blue",   linestyle=":",  alpha=0.7, label="LE=0.3 ○ 良好")
    ax.axvline(x=0.2, color="orange", linestyle="-.", alpha=0.7, label="LE=0.2 △ 許容")
    ax.legend(fontsize=9)

    ax.invert_yaxis()
    ax.grid(axis="x", alpha=0.3)

    for i, (bar, e) in enumerate(zip(bars, entries)):
        ax.text(e["le"] + 0.005, i, f"{e['le']:.4f}", va="center", fontsize=8)

    plt.tight_layout()
    out_path = out_dir / "adcp_vs_smina_le_comparison.png"
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  LE 比較グラフ保存: {out_path}")


# ──────────────────────────────────────────────────────────────────────────────
# レポート・サマリー
# ──────────────────────────────────────────────────────────────────────────────

def write_report(adcp: ADCPResult, smina: List[dict], out_dir: Path):
    """JSON レポートを書き出す"""
    report = {
        "adcp_cyclic_peptide": {
            "sequence":          adcp.sequence,
            "affinity_kcal_mol": adcp.affinity if adcp.affinity < float("inf") else None,
            "best_energy":       adcp.best_energy if adcp.best_energy < float("inf") else None,
            "n_clusters":        adcp.n_clusters,
            "hac":               adcp.hac,
            "le_adcp":           adcp.le,
            "status":            adcp.status,
            "dlg_file":          adcp.dlg_file,
            "pdb_file":          adcp.pdb_file,
        },
        "smina_top_small_molecules": [
            {
                "rank":         i + 1,
                "name":         r["name"],
                "score_kcal_mol": r["score"],
                "hac":          r.get("hac"),
                "le_smina":     r.get("le"),
                "pair":         r.get("pair", ""),
            }
            for i, r in enumerate(smina)
        ],
        "comparison_note": (
            "ADCP (AutoDock CrankPep) と smina は異なるスコアリング関数を使用しています。"
            "ADCP は Monte Carlo 法 + implicit 溶媒で環状ペプチドに特化しており、"
            "smina は AutoDock Vina スコア関数で低分子に使用します。"
            "スコアの直接比較はできません。LE も参考値として扱ってください。"
        ),
        "scoring_info": {
            "ADCP": {
                "tool":     "AutoDock CrankPep (adcp v0.0.25)",
                "method":   "Monte Carlo + OpenMM implicit solvent, cyclic mode (-cyc)",
                "typical_range": "-10 ~ -35 kcal/mol",
            },
            "smina": {
                "tool":     "smina (AutoDock Vina based)",
                "method":   "AutoDock Vina scoring function",
                "typical_range": "-5 ~ -10 kcal/mol",
            },
        },
    }

    out_path = out_dir / "adcp_comparison_report.json"
    with open(out_path, "w", encoding="utf-8") as f:
        json.dump(report, f, indent=2, ensure_ascii=False)
    print(f"  レポート保存: {out_path}")


def print_summary(adcp: ADCPResult, smina: List[dict]):
    """コンソールにサマリーを表示"""
    print()
    print("=" * 80)
    print("  ADCP CrankPep — 環状ペプチド ドッキング結果")
    print("=" * 80)
    if adcp.status == "success" and adcp.affinity < float("inf"):
        le_str = f"{adcp.le:.4f}" if adcp.le else "N/A"
        print(f"  配列         : cyclo-{adcp.sequence}")
        print(f"  Affinity     : {adcp.affinity:.3f} kcal/mol  (ADCP スコア)")
        print(f"  Best energy  : {adcp.best_energy:.3f} kcal/mol")
        print(f"  HAC          : {adcp.hac or 'N/A'}")
        print(f"  LE (ADCP)    : {le_str}")
        print(f"  Clusters     : {adcp.n_clusters}")
        print(f"  DLG          : {adcp.dlg_file}")
        print(f"  PDB          : {adcp.pdb_file or 'N/A'}")
    else:
        print(f"  ドッキング失敗: status={adcp.status}")
        print(f"  エラー: {adcp.error_message}")

    print()
    print("=" * 80)
    print(f"  smina — 低分子 ドッキング結果 Top {len(smina)}")
    print("=" * 80)
    print(f"  {'順位':<4}  {'スコア(smina)':>14}  {'HAC':>5}  {'LE(smina)':>10}  分子名")
    print("  " + "-" * 76)
    for i, r in enumerate(smina, 1):
        le_str = f"{r['le']:.4f}" if r.get("le") else "   N/A"
        print(f"  {i:<4}  {r['score']:>14.3f}  "
              f"{str(r.get('hac','?')):>5}  {le_str:>10}  {r['name'][:45]}")

    print()
    print("─" * 80)
    print("  ⚠  スコアスケールについて")
    print("     ADCP  : Monte Carlo + implicit 溶媒 (典型: -10〜-35 kcal/mol)")
    print("     smina : AutoDock Vina スコア      (典型: -5〜-10 kcal/mol)")
    print("     → 2手法のスコアを直接比較することはできません。")
    print("     → LE も参考値として扱ってください。")
    print("─" * 80)


# ──────────────────────────────────────────────────────────────────────────────
# メイン
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=(
            "AutoDock CrankPep で環状ペプチドをドッキングし、"
            "smina 低分子結果と比較する。"
        ),
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "--sequence", default=DEFAULT_SEQUENCE,
        help=f"ペプチド配列 (デフォルト: {DEFAULT_SEQUENCE})",
    )
    parser.add_argument(
        "--setup-receptor", action="store_true",
        help="受容体ターゲットファイル (.trg) を再作成する",
    )
    parser.add_argument(
        "--quick", action="store_true",
        help="クイックモード: N=1, n=5000 (テスト用、精度は低い)",
    )
    parser.add_argument(
        "--n-runs", type=int, default=5,
        help="MC 探索回数 (デフォルト: 5)",
    )
    parser.add_argument(
        "--n-evals", type=int, default=100000,
        help="評価ステップ数 (デフォルト: 100000)",
    )
    parser.add_argument(
        "--timeout", type=int, default=3600,
        help="ドッキングのタイムアウト秒数 (デフォルト: 3600)",
    )
    parser.add_argument(
        "--force", action="store_true",
        help="既存結果を上書きして再実行する",
    )
    args = parser.parse_args()

    # クイックモード上書き
    if args.quick:
        args.n_runs  = 1
        args.n_evals = 5000
        args.timeout = 180
        print("  [クイックモード] N=1, n=5000, timeout=180s (精度は低い)")

    print("\n" + "=" * 65)
    print("  AutoDock CrankPep ドッキング + smina 比較")
    print("=" * 65)

    # ① micromamba 確認
    mamba = find_micromamba()
    if mamba is None:
        print("[エラー] micromamba が見つかりません。")
        print("  インストール方法: curl -Ls https://micro.mamba.pm/install.sh | bash")
        sys.exit(1)
    print(f"  micromamba : {mamba}")

    # ② adcpsuite 環境確認
    if not check_adcpsuite_env(mamba):
        print("[エラー] adcpsuite 環境が見つかりません。")
        print("  セットアップ: bash adcpsuite_micromamba.sh")
        sys.exit(1)
    print("  adcpsuite  : OK")

    ADCP_DIR.mkdir(parents=True, exist_ok=True)

    # ③ ターゲットファイル準備
    print("\n  ① 受容体ターゲットファイル準備 ...")
    try:
        trg_file = setup_receptor(mamba, force=args.setup_receptor)
    except (FileNotFoundError, RuntimeError) as e:
        print(f"[エラー] {e}")
        sys.exit(1)

    # ④ ADCP ドッキング
    print(f"\n  ② ADCP ドッキング実行  (cyclo-{args.sequence})")
    print(f"     N={args.n_runs}, n={args.n_evals}, timeout={args.timeout}s")
    adcp_result = run_adcp_docking(
        mamba    = mamba,
        sequence = args.sequence,
        trg_file = trg_file,
        n_runs   = args.n_runs,
        n_evals  = args.n_evals,
        timeout  = args.timeout,
        force    = args.force,
    )
    print(f"     ステータス: {adcp_result.status}")
    if adcp_result.status == "success":
        print(f"     Affinity : {adcp_result.affinity:.3f} kcal/mol")
    else:
        print(f"     エラー   : {adcp_result.error_message}")

    # ⑤ smina 結果読み込み
    print("\n  ③ smina 低分子結果読み込み中 ...")
    smina_records = load_smina_results(top_n=15)
    print(f"     {len(smina_records)} 件ロード")

    # ⑥ 比較グラフ
    print("\n  ④ 比較グラフ生成中 ...")
    plot_comparison(adcp_result, smina_records, ADCP_DIR)
    plot_le_comparison(adcp_result, smina_records, ADCP_DIR)

    # ⑦ レポート保存
    write_report(adcp_result, smina_records, ADCP_DIR)

    # ⑧ サマリー表示
    print_summary(adcp_result, smina_records)

    print(f"\n  完了! 出力先: {ADCP_DIR}")


if __name__ == "__main__":
    main()
