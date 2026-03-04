"""
analyze_prodigy.py
==================
PRODIGY (PROtein binDIng enerGY prediction) アルゴリズムの手実装。

論文: Vangone & Bonvin, eLife (2015)
      Xue et al., Bioinformatics (2016)

アルゴリズム:
  1. 界面接触を残基タイプ (C/P/N) に分類
     C = Charged  (ASP, GLU, LYS, ARG, HIS)
     P = Polar    (SER, THR, ASN, GLN, TYR, TRP, CYS)
     N = Non-polar(ALA, VAL, ILE, LEU, MET, PHE, PRO, GLY)

  2. 界面接触を6タイプにカウント (Cβ間距離 ≤ 5.5 Å, GLYはCα)
     CC, CP, CN, PP, PN, NN

  3. NIS (Non-Interface Surface) 残基をタイプ別にカウント
     → 溶媒に露出した非界面残基の組成が親和性に影響

  4. 線形回帰式で ΔG を予測 (単位: kcal/mol)
     ΔG = w_CC×IC_CC + w_CP×IC_CP + w_CN×IC_CN
        + w_PP×IC_PP + w_PN×IC_PN + w_NN×IC_NN
        + w_NISP×NIS_P + w_NISC×NIS_C + intercept

  5. Kd を計算: Kd = exp(ΔG / RT)

参照係数: Vangone & Bonvin (2015) Table 1
"""

import math
import numpy as np
from Bio.PDB import PDBParser, NeighborSearch
import warnings
warnings.filterwarnings("ignore")


# ──────────────────────────────────────────
# 残基タイプ分類
# ──────────────────────────────────────────
CHARGED_AA   = {"ASP", "GLU", "LYS", "ARG", "HIS"}
POLAR_AA     = {"SER", "THR", "ASN", "GLN", "TYR", "TRP", "CYS"}
NONPOLAR_AA  = {"ALA", "VAL", "ILE", "LEU", "MET", "PHE", "PRO", "GLY"}

def residue_type(resname: str) -> str:
    """残基をC/P/Nに分類"""
    if resname in CHARGED_AA:   return "C"
    if resname in POLAR_AA:     return "P"
    return "N"  # non-polar or unknown


# ──────────────────────────────────────────
# PRODIGY 回帰係数 (Vangone & Bonvin 2015)
# ──────────────────────────────────────────
PRODIGY_WEIGHTS = {
    "IC_CC" : -0.09459,
    "IC_CP" : -0.10007,
    "IC_CN" :  0.19577,   # charged - nonpolar
    "IC_PP" : -0.22671,
    "IC_PN" :  0.18681,   # polar - nonpolar
    "IC_NN" :  0.34282,
    "NIS_P" :  0.01384,   # 非界面極性残基
    "NIS_C" :  0.02064,   # 非界面荷電残基
    "intercept": -15.9433,
}

TEMPERATURE  = 298.15   # K
GAS_CONSTANT = 1.987e-3  # kcal / (mol·K)


# ──────────────────────────────────────────
# 界面接触のカウント
# ──────────────────────────────────────────

def _get_cb_atom(residue):
    """Cβ原子を返す（GLYはCα）"""
    if "CB" in [a.get_name() for a in residue.get_atoms()]:
        return residue["CB"]
    if "CA" in [a.get_name() for a in residue.get_atoms()]:
        return residue["CA"]
    return None


def count_interface_contacts(model,
                              chain_a: str,
                              chain_b: str,
                              cutoff: float = 5.5) -> dict:
    """
    Cβ/Cα 間距離に基づく界面接触を6タイプでカウント。

    Returns:
        {"IC_CC": n, "IC_CP": n, "IC_CN": n,
         "IC_PP": n, "IC_PN": n, "IC_NN": n,
         "contacts": [(res_a, type_a, res_b, type_b, dist), ...]}
    """
    residues_a = list(model[chain_a].get_residues())
    residues_b = list(model[chain_b].get_residues())

    counts = {"IC_CC": 0, "IC_CP": 0, "IC_CN": 0,
              "IC_PP": 0, "IC_PN": 0, "IC_NN": 0}
    contacts = []

    for res_a in residues_a:
        cb_a = _get_cb_atom(res_a)
        if cb_a is None:
            continue
        ta = residue_type(res_a.get_resname())

        for res_b in residues_b:
            cb_b = _get_cb_atom(res_b)
            if cb_b is None:
                continue
            tb = residue_type(res_b.get_resname())

            dist = float(np.linalg.norm(cb_a.coord - cb_b.coord))
            if dist > cutoff:
                continue

            # タイプペアを正規化 (アルファベット順)
            pair = "".join(sorted([ta, tb]))
            key  = f"IC_{pair}"
            if key in counts:
                counts[key] += 1

            contacts.append({
                "res_a"  : res_a.get_resname(),
                "resnum_a": res_a.get_id()[1],
                "type_a" : ta,
                "res_b"  : res_b.get_resname(),
                "resnum_b": res_b.get_id()[1],
                "type_b" : tb,
                "distance": round(dist, 3),
            })

    return {**counts, "contacts": contacts}


# ──────────────────────────────────────────
# NIS (Non-Interface Surface) 残基のカウント
# ──────────────────────────────────────────

def count_nis_residues(model,
                        chain_a: str,
                        chain_b: str,
                        interface_contacts: dict,
                        sasa_result: dict | None = None) -> dict:
    """
    NIS残基 = 表面に露出しているが界面に参加していない残基。

    sasa_resultがある場合はSASA > 0 を「表面露出」の基準とする。
    ない場合はすべての残基をNIS候補とみなす（近似）。

    Returns: {"NIS_P": n, "NIS_C": n}
    """
    # 界面に参加している残基セット
    iface_a = set(c["resnum_a"] for c in interface_contacts["contacts"])
    iface_b = set(c["resnum_b"] for c in interface_contacts["contacts"])

    nis_p = 0
    nis_c = 0

    for chain_id, iface_set in [(chain_a, iface_a), (chain_b, iface_b)]:
        for residue in model[chain_id].get_residues():
            resnum  = residue.get_id()[1]
            resname = residue.get_resname()
            if resnum in iface_set:
                continue   # 界面残基は除外

            # 表面露出チェック
            if sasa_result:
                key = (chain_id, resnum, resname)
                alone_key = "protein_sasa" if chain_id == chain_a else "peptide_sasa"
                sasa_val  = sasa_result[alone_key].get(key, 0)
                if sasa_val < 1.0:
                    continue  # 埋没残基はNIS対象外

            rtype = residue_type(resname)
            if rtype == "P":
                nis_p += 1
            elif rtype == "C":
                nis_c += 1

    return {"NIS_P": nis_p, "NIS_C": nis_c}


# ──────────────────────────────────────────
# ΔG・Kd 計算
# ──────────────────────────────────────────

def predict_binding_affinity(ic_counts: dict,
                              nis_counts: dict,
                              temperature: float = TEMPERATURE) -> dict:
    """
    PRODIGY 線形回帰式でΔGとKdを予測する。

    Args:
        ic_counts : {"IC_CC": n, "IC_CP": n, ...}
        nis_counts: {"NIS_P": n, "NIS_C": n}
        temperature: K
    Returns:
        {"dG": kcal/mol, "Kd": M, "Kd_str": "xnM/μM/mM"}
    """
    w = PRODIGY_WEIGHTS
    dG = (w["IC_CC"] * ic_counts.get("IC_CC", 0) +
          w["IC_CP"] * ic_counts.get("IC_CP", 0) +
          w["IC_CN"] * ic_counts.get("IC_CN", 0) +
          w["IC_PP"] * ic_counts.get("IC_PP", 0) +
          w["IC_PN"] * ic_counts.get("IC_PN", 0) +
          w["IC_NN"] * ic_counts.get("IC_NN", 0) +
          w["NIS_P"] * nis_counts.get("NIS_P", 0) +
          w["NIS_C"] * nis_counts.get("NIS_C", 0) +
          w["intercept"])

    kd = math.exp(dG / (GAS_CONSTANT * temperature))

    # Kd を読みやすい単位に変換
    if kd < 1e-9:
        kd_str = f"{kd * 1e12:.1f} pM"
    elif kd < 1e-6:
        kd_str = f"{kd * 1e9:.1f} nM"
    elif kd < 1e-3:
        kd_str = f"{kd * 1e6:.1f} μM"
    else:
        kd_str = f"{kd * 1e3:.1f} mM"

    return {"dG": round(dG, 3), "Kd": kd, "Kd_str": kd_str}


# ──────────────────────────────────────────
# メイン解析関数
# ──────────────────────────────────────────

def run_prodigy(pdb_path: str,
                chain_protein: str = "A",
                chain_peptide: str = "B",
                sasa_result: dict | None = None,
                cutoff: float = 5.5) -> dict:
    """
    PRODIGY 解析を実行して ΔG・Kd を予測する。

    Args:
        pdb_path     : 複合体PDBファイルパス
        chain_protein: タンパク質チェーンID
        chain_peptide: ペプチドチェーンID
        sasa_result  : analyze_sasa.calc_delta_sasa() の結果 (オプション)
        cutoff       : Cβ/Cα 接触距離カットオフ Å (default 5.5)
    Returns:
        dict with full results
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("complex", pdb_path)
    model = structure[0]

    ic  = count_interface_contacts(model, chain_protein, chain_peptide, cutoff)
    nis = count_nis_residues(model, chain_protein, chain_peptide, ic, sasa_result)
    aff = predict_binding_affinity(ic, nis)

    return {
        "interface_contacts": {k: v for k, v in ic.items() if k != "contacts"},
        "contact_details"   : ic["contacts"],
        "nis_residues"      : nis,
        "dG_kcal_mol"       : aff["dG"],
        "Kd"                : aff["Kd"],
        "Kd_str"            : aff["Kd_str"],
        "total_ic"          : sum(v for k, v in ic.items() if k.startswith("IC_")),
    }


def print_prodigy_summary(result: dict):
    print("\n" + "=" * 65)
    print("  PRODIGY 結合親和性予測")
    print("=" * 65)
    print(f"  予測 ΔG    : {result['dG_kcal_mol']:.3f} kcal/mol")
    print(f"  予測 Kd    : {result['Kd_str']}")
    print(f"  総界面接触数: {result['total_ic']}")
    print()
    print("  界面接触タイプ内訳:")
    ic = result["interface_contacts"]
    type_names = {
        "IC_CC": "Charged–Charged (CC)",
        "IC_CP": "Charged–Polar   (CP)",
        "IC_CN": "Charged–NonPolar(CN)",
        "IC_PP": "Polar–Polar     (PP)",
        "IC_PN": "Polar–NonPolar  (PN)",
        "IC_NN": "NonPolar–NonPolar(NN)",
    }
    for key, label in type_names.items():
        n   = ic.get(key, 0)
        bar = "█" * n
        print(f"  {label}: {n:>3}  {bar}")
    print()
    print(f"  NIS残基 (非界面表面残基):")
    print(f"    極性   (P): {result['nis_residues']['NIS_P']}")
    print(f"    荷電   (C): {result['nis_residues']['NIS_C']}")
    print("=" * 65)


def save_prodigy_json(result: dict, output_path: str):
    import json
    out = {k: v for k, v in result.items() if k != "contact_details"}
    out["Kd"] = float(result["Kd"])
    with open(output_path, "w") as f:
        json.dump(out, f, indent=2, ensure_ascii=False)
    print(f"  PRODIGY 結果 JSON 保存: {output_path}")


def plot_prodigy_contacts(result: dict, output_path: str):
    """界面接触タイプの円グラフ"""
    import plot_utils  # noqa: F401
    import matplotlib.pyplot as plt

    ic     = result["interface_contacts"]
    labels = ["CC", "CP", "CN", "PP", "PN", "NN"]
    keys   = [f"IC_{l}" for l in labels]
    values = [ic.get(k, 0) for k in keys]
    colors = ["#e74c3c", "#e67e22", "#f1c40f",
              "#2ecc71", "#3498db", "#9b59b6"]

    # ゼロ値を除外
    filtered = [(l, v, c) for l, v, c in zip(labels, values, colors) if v > 0]
    if not filtered:
        return
    labels_f  = [f[0] for f in filtered]
    values_f  = [f[1] for f in filtered]
    colors_f  = [f[2] for f in filtered]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # 円グラフ
    wedges, texts, autotexts = ax1.pie(
        values_f, labels=labels_f, colors=colors_f,
        autopct="%1.0f%%", startangle=90,
        textprops={"fontsize": 11}
    )
    ax1.set_title("界面接触タイプ分布", fontsize=12, fontweight="bold")

    # 棒グラフ + 係数
    weights = [abs(PRODIGY_WEIGHTS[f"IC_{l}"]) for l in labels_f]
    contributions = [v * w for v, w in zip(values_f, weights)]
    bars = ax2.bar(labels_f, contributions, color=colors_f, edgecolor="white")
    ax2.set_xlabel("接触タイプ", fontsize=11)
    ax2.set_ylabel("|係数 × 接触数|", fontsize=11)
    ax2.set_title("ΔG への寄与 (絶対値)", fontsize=12, fontweight="bold")

    dg   = result["dG_kcal_mol"]
    kdst = result["Kd_str"]
    fig.suptitle(f"PRODIGY: ΔG = {dg:.2f} kcal/mol  /  Kd = {kdst}",
                 fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  PRODIGY グラフ 保存: {output_path}")


if __name__ == "__main__":
    import sys
    pdb = sys.argv[1] if len(sys.argv) > 1 else "Protein_Peptide.pdb"
    result = run_prodigy(pdb)
    print_prodigy_summary(result)
    save_prodigy_json(result, "prodigy_result.json")
    plot_prodigy_contacts(result, "prodigy_contacts.png")
