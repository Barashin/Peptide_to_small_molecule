"""
visualize.py
============
相互作用解析・ファーマコフォア・候補分子の可視化。
"""

import numpy as np
import plot_utils  # noqa: F401  日本語フォント設定
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import Counter, defaultdict


FEAT_COLOR = {
    "HBA": "#e74c3c",   # 赤
    "HBD": "#3498db",   # 青
    "HYD": "#f1c40f",   # 黄
    "ARO": "#e67e22",   # オレンジ
    "POS": "#1abc9c",   # シアン
    "NEG": "#9b59b6",   # 紫
}

CONTACT_COLOR = {
    "hbond"        : "#e74c3c",
    "electrostatic": "#9b59b6",
    "hydrophobic"  : "#f1c40f",
    "van_der_waals": "#95a5a6",
    "other"        : "#bdc3c7",
}


def plot_residue_scores(residue_scores: dict, output_path: str):
    """ペプチド残基の重要度スコアを棒グラフで可視化"""
    if not residue_scores:
        return

    labels = [f"{name}{num}" for (num, name) in residue_scores.keys()]
    values = list(residue_scores.values())

    fig, ax = plt.subplots(figsize=(max(8, len(labels) * 0.9), 5))
    bars = ax.bar(labels, values, color="#2ecc71", edgecolor="white", linewidth=0.8)

    # スコアが高い上位3残基を強調
    sorted_vals = sorted(enumerate(values), key=lambda x: -x[1])
    for rank, (i, v) in enumerate(sorted_vals[:3]):
        bars[i].set_color(["#e74c3c", "#e67e22", "#f1c40f"][rank])

    ax.set_xlabel("ペプチド残基", fontsize=12)
    ax.set_ylabel("相互作用スコア", fontsize=12)
    ax.set_title("タンパク質-ペプチド 残基別相互作用スコア", fontsize=14, fontweight="bold")
    ax.set_ylim(0, max(values) * 1.2)

    # 凡例
    patches = [mpatches.Patch(color=c, label=l) for c, l in
               [("#e74c3c","1位"), ("#e67e22","2位"), ("#f1c40f","3位"), ("#2ecc71","その他")]]
    ax.legend(handles=patches, loc="upper right", framealpha=0.8)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  残基スコア図 保存: {output_path}")


def plot_interaction_map(contacts: list[dict], output_path: str):
    """タンパク質残基 × ペプチド残基 の相互作用マップ"""
    from analyze_interactions import classify_contact

    # 接触タイプ別にカウント
    matrix = defaultdict(Counter)
    for c in contacts:
        ct  = classify_contact(c)
        matrix[(c["protein_resnum"], c["protein_res"])][(c["peptide_resnum"], c["peptide_res"])] += 1

    prot_keys = sorted(set(c["protein_resnum"] for c in contacts))
    pep_keys  = sorted(set(c["peptide_resnum"] for c in contacts))

    if not prot_keys or not pep_keys:
        return

    # ヒートマップ行列
    data = np.zeros((len(prot_keys), len(pep_keys)))
    for c in contacts:
        pi = prot_keys.index(c["protein_resnum"])
        pj = pep_keys.index(c["peptide_resnum"])
        data[pi, pj] += 1

    # 接触が多い上位20残基に絞る (大きな構造の場合)
    row_sums = data.sum(axis=1)
    top_rows = np.argsort(row_sums)[-20:][::-1]
    data_trimmed = data[top_rows, :]
    prot_labels_trimmed = [f"{prot_keys[i]}" for i in top_rows]

    # 対応する残基名も取得
    prot_res_map = {c["protein_resnum"]: c["protein_res"] for c in contacts}
    prot_labels_trimmed = [f"{prot_res_map.get(prot_keys[i], '')}{prot_keys[i]}"
                           for i in top_rows]
    pep_res_map   = {c["peptide_resnum"]: c["peptide_res"] for c in contacts}
    pep_labels    = [f"{pep_res_map.get(k, '')}{k}" for k in pep_keys]

    fig, ax = plt.subplots(figsize=(max(8, len(pep_keys) * 1.2),
                                    max(6, len(top_rows) * 0.5)))
    im = ax.imshow(data_trimmed, cmap="YlOrRd", aspect="auto")
    ax.set_xticks(range(len(pep_labels)))
    ax.set_xticklabels(pep_labels, rotation=45, ha="right", fontsize=10)
    ax.set_yticks(range(len(prot_labels_trimmed)))
    ax.set_yticklabels(prot_labels_trimmed, fontsize=8)
    ax.set_xlabel("ペプチド残基", fontsize=12)
    ax.set_ylabel("タンパク質残基 (上位20)", fontsize=12)
    ax.set_title("相互作用マップ (接触数)", fontsize=14, fontweight="bold")
    plt.colorbar(im, ax=ax, shrink=0.8, label="接触数")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  相互作用マップ 保存: {output_path}")


def plot_pharmacophore_3d(features: list[dict], output_path: str):
    """ファーマコフォア特徴点の3D散布図"""
    if not features:
        return

    fig = plt.figure(figsize=(10, 8))
    ax  = fig.add_subplot(111, projection="3d")

    grouped = defaultdict(list)
    for f in features:
        grouped[f["type"]].append(f["coord"])

    handles = []
    for feat_type, coords in grouped.items():
        coords = np.array(coords)
        color  = FEAT_COLOR.get(feat_type, "gray")
        sc = ax.scatter(coords[:, 0], coords[:, 1], coords[:, 2],
                        c=color, s=120, alpha=0.85, label=feat_type, edgecolors="white")
        handles.append(sc)

    # 残基ラベル
    for f in features:
        x, y, z = f["coord"]
        ax.text(x, y, z + 0.3, f"{f['resname']}{f['resnum']}", fontsize=6, alpha=0.7)

    ax.set_xlabel("X (Å)")
    ax.set_ylabel("Y (Å)")
    ax.set_zlabel("Z (Å)")
    ax.set_title("ファーマコフォアモデル (3D)", fontsize=13, fontweight="bold")
    ax.legend(loc="upper left", fontsize=9)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    print(f"  ファーマコフォア3D図 保存: {output_path}")


def plot_candidate_structures(candidates: list[dict], output_path: str):
    """RDKit で候補分子の2D構造画像を生成"""
    from rdkit.Chem import Draw
    from rdkit.Chem import AllChem

    mols   = [c["mol"] for c in candidates if c["mol"] is not None]
    names  = [c["name"][:30] for c in candidates if c["mol"] is not None]

    if not mols:
        return

    n = len(mols)
    cols = min(n, 3)
    rows = (n + cols - 1) // cols

    img = Draw.MolsToGridImage(
        mols, molsPerRow=cols, subImgSize=(400, 300),
        legends=names, returnPNG=False
    )
    img.save(output_path)
    print(f"  候補分子構造図 保存: {output_path}")
