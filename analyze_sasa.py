"""
analyze_sasa.py
===============
Solvent Accessible Surface Area (SASA) 差分解析モジュール。

複合体・タンパク質単体・ペプチド単体それぞれのSASAを計算し、
ΔSASAが大きい残基（結合によって埋没した残基）をホットスポットとして同定する。

アルゴリズム:
  1. 複合体全体のSASAを計算
  2. タンパク質チェーン単体のSASAを計算
  3. ペプチドチェーン単体のSASAを計算
  4. ΔSASA = SASA(単体) - SASA(複合体) を各残基で算出
  5. ΔSASAが大きい残基を界面残基・ホットスポットとして報告

参照: Lee & Richards (1971), Shrake & Rupley (1973)
"""

import numpy as np
from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
import io
import warnings
warnings.filterwarnings("ignore")


# ΔSASAのホットスポット閾値 (Å²)
HOTSPOT_THRESHOLD = 1.0    # 界面残基の最小ΔSASA
CORE_THRESHOLD    = 10.0   # コア界面残基の閾値


class ChainSelect(Select):
    """指定チェーンのみを選択するセレクタ"""
    def __init__(self, chain_ids):
        self.chain_ids = set(chain_ids)

    def accept_chain(self, chain):
        return chain.id in self.chain_ids


def _structure_to_handle(structure, chain_ids: list) -> io.StringIO:
    """構造を指定チェーンのみPDB文字列として返す"""
    pdbio = PDBIO()
    pdbio.set_structure(structure)
    handle = io.StringIO()
    pdbio.save(handle, ChainSelect(chain_ids))
    handle.seek(0)
    return handle


def _load_from_handle(handle, name: str):
    """StringIOからBioPython構造を読み込む"""
    parser = PDBParser(QUIET=True)
    return parser.get_structure(name, handle)


def compute_sasa(structure, level: str = "R", probe_radius: float = 1.4,
                 n_points: int = 100) -> dict:
    """
    ShrakeRupley アルゴリズムでSASAを計算。

    Args:
        structure : BioPython structure object
        level     : "R" (残基) or "A" (原子)
        probe_radius: 探針半径 Å (水=1.4)
        n_points  : 球面サンプリング点数 (多いほど精度↑)
    Returns:
        dict: {(chain_id, res_id): sasa_value}
    """
    sr = ShrakeRupley(probe_radius=probe_radius, n_points=n_points)
    sr.compute(structure, level=level)

    result = {}
    for model in structure:
        for chain in model:
            for residue in chain:
                key = (chain.id, residue.get_id()[1], residue.get_resname())
                result[key] = residue.sasa
    return result


def calc_delta_sasa(pdb_path: str,
                    chain_protein: str = "A",
                    chain_peptide: str = "B",
                    probe_radius: float = 1.4,
                    n_points: int = 100) -> dict:
    """
    複合体と各単体のSASAを比較してΔSASAを計算する。

    Returns:
        dict with keys:
          "complex_sasa"  : {(chain, resnum, resname): sasa}  複合体
          "protein_sasa"  : {(chain, resnum, resname): sasa}  タンパク質単体
          "peptide_sasa"  : {(chain, resnum, resname): sasa}  ペプチド単体
          "delta_protein" : {(chain, resnum, resname): Δsasa}  タンパク質側の埋没量
          "delta_peptide" : {(chain, resnum, resname): Δsasa}  ペプチド側の埋没量
          "interface_residues_protein": list of (resnum, resname, delta)
          "interface_residues_peptide": list of (resnum, resname, delta)
          "total_buried_sasa": float  全埋没面積の合計 Å²
    """
    parser = PDBParser(QUIET=True)
    struct_complex = parser.get_structure("complex", pdb_path)

    # 各単体構造を作成
    handle_prot = _structure_to_handle(struct_complex, [chain_protein])
    handle_pep  = _structure_to_handle(struct_complex, [chain_peptide])
    struct_prot = _load_from_handle(handle_prot, "protein")
    struct_pep  = _load_from_handle(handle_pep,  "peptide")

    # SASA計算（複合体・単体それぞれ）
    sasa_complex = compute_sasa(struct_complex, probe_radius=probe_radius, n_points=n_points)
    sasa_protein = compute_sasa(struct_prot,    probe_radius=probe_radius, n_points=n_points)
    sasa_peptide = compute_sasa(struct_pep,     probe_radius=probe_radius, n_points=n_points)

    # ΔSASAの計算 (単体 - 複合体)
    delta_protein = {}
    for key, sasa_alone in sasa_protein.items():
        sasa_bound = sasa_complex.get(key, sasa_alone)
        delta_protein[key] = sasa_alone - sasa_bound  # 正の値 = 複合体で埋没

    delta_peptide = {}
    for key, sasa_alone in sasa_peptide.items():
        sasa_bound = sasa_complex.get(key, sasa_alone)
        delta_peptide[key] = sasa_alone - sasa_bound

    # 界面残基の同定
    iface_protein = sorted(
        [(key[1], key[2], d) for key, d in delta_protein.items() if d >= HOTSPOT_THRESHOLD],
        key=lambda x: -x[2]
    )
    iface_peptide = sorted(
        [(key[1], key[2], d) for key, d in delta_peptide.items() if d >= HOTSPOT_THRESHOLD],
        key=lambda x: -x[2]
    )

    total_buried = sum(d for d in delta_protein.values() if d > 0) + \
                   sum(d for d in delta_peptide.values() if d > 0)

    return {
        "complex_sasa"               : sasa_complex,
        "protein_sasa"               : sasa_protein,
        "peptide_sasa"               : sasa_peptide,
        "delta_protein"              : delta_protein,
        "delta_peptide"              : delta_peptide,
        "interface_residues_protein" : iface_protein,
        "interface_residues_peptide" : iface_peptide,
        "total_buried_sasa"          : round(total_buried, 2),
    }


def score_residues_by_sasa(delta_peptide: dict) -> dict:
    """
    ΔSASAベースの残基スコアを返す。
    distance-baseのスコアと統合するために {(resnum, resname): score} 形式で返す。
    """
    scores = {}
    for (chain, resnum, resname), delta in delta_peptide.items():
        if delta > 0:
            scores[(resnum, resname)] = round(delta, 3)
    return dict(sorted(scores.items(), key=lambda x: -x[1]))


def print_sasa_summary(result: dict):
    print("\n" + "=" * 65)
    print("  ΔSASA 解析結果")
    print("=" * 65)
    print(f"  総埋没面積 (BSA): {result['total_buried_sasa']:.1f} Å²")
    print(f"  ※ タンパク質-タンパク質複合体の典型値: 1500〜3000 Å²")
    print()

    print("  ペプチド側 界面残基 (ΔSASA ≥ 1.0 Å²):")
    print(f"  {'残基':<10} {'ΔSASA (Å²)':>12}  {'役割'}")
    print("  " + "-" * 50)
    for resnum, resname, delta in result["interface_residues_peptide"]:
        role = "コア界面" if delta >= CORE_THRESHOLD else "周辺界面"
        bar  = "█" * int(delta / 5)
        print(f"  {resname}{resnum:<6} {delta:>12.1f}  {role}  {bar}")

    print()
    print(f"  タンパク質側 界面残基 (上位10):")
    print(f"  {'残基':<10} {'ΔSASA (Å²)':>12}")
    print("  " + "-" * 30)
    for resnum, resname, delta in result["interface_residues_protein"][:10]:
        print(f"  {resname}{resnum:<6} {delta:>12.1f}")
    print("=" * 65)


def save_sasa_csv(result: dict, output_path: str):
    import csv
    with open(output_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["chain", "resnum", "resname",
                         "sasa_alone", "sasa_complex", "delta_sasa", "is_interface"])
        # ペプチド残基
        for (chain, resnum, resname), delta in result["delta_peptide"].items():
            alone   = result["peptide_sasa"].get((chain, resnum, resname), 0)
            complex_ = alone - delta
            writer.writerow([chain, resnum, resname,
                             round(alone, 3), round(complex_, 3),
                             round(delta, 3), delta >= HOTSPOT_THRESHOLD])
        # タンパク質残基
        for (chain, resnum, resname), delta in result["delta_protein"].items():
            alone   = result["protein_sasa"].get((chain, resnum, resname), 0)
            complex_ = alone - delta
            writer.writerow([chain, resnum, resname,
                             round(alone, 3), round(complex_, 3),
                             round(delta, 3), delta >= HOTSPOT_THRESHOLD])
    print(f"  ΔSASA CSV 保存: {output_path}")


def plot_sasa_comparison(result: dict, output_path: str):
    """ペプチド残基のSASA比較（単体 vs 複合体）棒グラフ"""
    import plot_utils  # noqa: F401
    import matplotlib.pyplot as plt

    iface = result["interface_residues_peptide"]
    if not iface:
        return

    labels = [f"{name}{num}" for num, name, _ in iface]
    deltas = [d for _, _, d in iface]
    colors = ["#e74c3c" if d >= CORE_THRESHOLD else "#e67e22" for d in deltas]

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.9), 5))
    ax.bar(labels, deltas, color=colors, edgecolor="white")
    ax.axhline(HOTSPOT_THRESHOLD, color="gray",  linestyle="--",
               linewidth=1, label=f"界面閾値 {HOTSPOT_THRESHOLD} Å²")
    ax.axhline(CORE_THRESHOLD,    color="#c0392b", linestyle="--",
               linewidth=1, label=f"コア界面閾値 {CORE_THRESHOLD} Å²")
    ax.set_xlabel("ペプチド残基", fontsize=12)
    ax.set_ylabel("ΔSASA (Å²)", fontsize=12)
    ax.set_title("結合界面の埋没表面積 (ΔSASA)", fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  ΔSASA グラフ 保存: {output_path}")


if __name__ == "__main__":
    import sys
    pdb = sys.argv[1] if len(sys.argv) > 1 else "Protein_Peptide.pdb"
    result = calc_delta_sasa(pdb)
    print_sasa_summary(result)
    save_sasa_csv(result, "sasa_analysis.csv")
    plot_sasa_comparison(result, "sasa_comparison.png")
