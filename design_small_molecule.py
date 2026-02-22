"""
design_small_molecule.py
========================
ペプチドの主要側鎖フラグメントを繋いで低分子リガンドを設計するモジュール。

アプローチ:
  1. 重要ペプチド残基の側鎖をRDKitフラグメントとしてマッピング
  2. 最も重要な2～4残基を選択
  3. ピペラジン/ベンゼン/ピリジンなどの骨格で連結
  4. 3D配座をファーマコフォア拘束つきで生成
  5. SDFとして出力

リンカー長の選択:
  - Cβ–Cβ距離（cbeta_coords を渡した場合）から自動計算
  - 座標がない場合は 短(C) / 中(CCCC) / 長(CCCCCC) の3パターンを生成
  - フラグメント自身の「リーチ」(Cβから官能基端までの推定距離)を考慮して
    実際に必要なリンカー長を補正する
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, rdMolDescriptors
from rdkit.Chem.rdchem import RWMol
from typing import Optional
import os

import linker_library as _ll

# FEgrow リンカー DB の遅延ロード
_LINKER_DB: list | None = None

def _get_linker_db() -> list:
    global _LINKER_DB
    if _LINKER_DB is None:
        _LINKER_DB = _ll.load_linker_db()
    return _LINKER_DB

def _get_linker_core(smiles_r1r2: str) -> str | None:
    """
    [R1]xxx[R2] 形式から線形コア SMILES を抽出する。
    [R1]CCC[R2] → 'CCC', [R1][R2] → ''
    環を含む非線形リンカーは None を返す。
    """
    s = smiles_r1r2
    # 先頭の [R1]/[R2] を除去
    if s.startswith("[R1]"):
        s = s[4:]
    elif s.startswith("[R2]"):
        s = s[4:]
    else:
        return None
    # 末尾の [R1]/[R2] を除去
    if s.endswith("[R2]"):
        s = s[:-4]
    elif s.endswith("[R1]"):
        s = s[:-4]
    else:
        return None
    # 残りに [R1]/[R2] がまだあれば非線形
    if "[R1]" in s or "[R2]" in s:
        return None
    return s  # "" は直接連結 (valid)


# 残基 → SMILES フラグメント (utils.residue_defs に統合)
from utils.residue_defs import RESIDUE_SMILES

# ──────────────────────────────────────────
# リンカー選択は FEgrow ライブラリ (linker_library.py) に委譲
# 座標なし時のフォールバック / 短距離時のフォールバック ([R1]/[R2] エントリ形式)
# ──────────────────────────────────────────
_FALLBACK_LINKER_ENTRIES = [
    {"name": "direct", "smiles_r1r2": "[R1][R2]",      "dist": 0.0, "hba": 0.0},
    {"name": "C4",     "smiles_r1r2": "[R1]CCCC[R2]",  "dist": 5.0, "hba": 0.0},
    {"name": "C6",     "smiles_r1r2": "[R1]CCCCCC[R2]","dist": 7.6, "hba": 0.0},
]
_SHORT_LINKER_ENTRIES = [
    {"name": "direct", "smiles_r1r2": "[R1][R2]",   "dist": 0.0, "hba": 0.0},
    {"name": "C1",     "smiles_r1r2": "[R1]C[R2]",  "dist": 1.3, "hba": 0.0},
    {"name": "C2",     "smiles_r1r2": "[R1]CC[R2]", "dist": 2.5, "hba": 0.0},
]
# 後方互換: strategy③ で使うコアSMILS (リスト)
_FALLBACK_LINKER_CORES = ["", "CCCC", "CCCCCC"]

from utils.residue_defs import FRAG_REACH

def select_linker(dist: float, reach_a: float, reach_b: float,
                  tolerance: float = 1.5) -> list[str]:
    """
    Cβ–Cβ距離とフラグメント両端のリーチから最適リンカーコアSMILESを選択。

    FEgrow ライブラリから |linker_dist - needed| <= tolerance のものを返す。
    必要リンカー長 = dist - reach_a - reach_b
    FEgrow 最小距離 (≈2.3 Å) 以下の場合は直接連結/短リンカーを返す。
    """
    needed = max(dist - reach_a - reach_b, 0.0)

    # FEgrow 最小距離より短い場合: シンプルな短リンカーで対応
    if needed <= 2.0:
        return ["", "C", "CC"][:max(1, int(needed / 0.8) + 1)]

    entries = _ll.select_linkers_by_dist(
        _get_linker_db(), needed, tolerance=tolerance, max_n=5,
    )
    cores: list[str] = []
    for e in entries:
        core = _get_linker_core(e["smiles_r1r2"])
        if core is not None and core not in cores:
            cores.append(core)
    if not cores:
        cores = ["", "CC", "CCCC"]  # 最終フォールバック
    return cores


# ──────────────────────────────────────────
# 低分子設計の骨格テンプレート
# key: 連結フラグメント数, value: (scaffold_name, SMARTS/SMILES with attachment points)
SCAFFOLD_TEMPLATES = {
    2: [
        ("piperazine",      "C1CN(CC([*:1])[*:2])CC1"),
        ("ethylene",        "[*:1]CC[*:2]"),
        ("phenyl_linker",   "[*:1]Cc1ccc(C[*:2])cc1"),
    ],
    3: [
        ("benzene_tri",     "c1c([*:1])cc([*:2])cc1[*:3]"),
        ("piperidine_tri",  "C1CC([*:1])CN(C1)[*:2].[*:3]"),   # dummy
        ("glycine_tri",     "[*:1]CC(=O)NCC([*:2])C(=O)NCC[*:3]"),
    ],
    4: [
        ("scaffold_4",      "[*:1]CC(=O)NCC([*:2])C(=O)NCC([*:3])C(=O)NCC[*:4]"),
    ],
}


def smiles_for_residue(res_name: str) -> Optional[str]:
    return RESIDUE_SMILES.get(res_name)


def select_key_residues(residue_scores: dict, top_n: int = 3) -> list[tuple]:
    """スコア上位 top_n 残基を選択 (GLYは除外)"""
    selected = []
    for (num, name), score in residue_scores.items():
        if name == "GLY":
            continue
        if smiles_for_residue(name) is None:
            continue
        selected.append((num, name, score))
        if len(selected) >= top_n:
            break
    return selected


def build_peptidomimetic(
    selected_residues: list[tuple],
    cbeta_coords: Optional[dict] = None,
) -> list[dict]:
    """
    選択した残基からペプチドミメティクスSMILESを複数生成。

    Args:
        selected_residues : [(resnum, resname, score), ...]  スコア降順
        cbeta_coords      : {resnum: np.ndarray(3)} — Cβ座標 (オプション)
                            渡すと残基間距離に合わせてリンカー長を自動選択。
                            None の場合は 短/中/長 の3パターンを生成。

    Returns:
        list of {'name', 'smiles', 'mol', 'source', 'linker'}
    """
    if len(selected_residues) == 0:
        return []

    # 各残基の側鎖SMILESを取得
    frags = []  # [(resnum, resname, clean_smiles)]
    for num, name, score in selected_residues:
        smi = smiles_for_residue(name)
        if smi:
            frags.append((num, name, smi.replace("[*]", "")))

    n_frags = len(frags)
    candidates = []
    seen_smiles = set()

    def _try_add(name, smiles_str, source, linker_name=""):
        """SMILESが有効なら重複除外してcandidatesに追加"""
        mol = Chem.MolFromSmiles(smiles_str)
        if mol is None:
            return
        canon = Chem.MolToSmiles(mol)
        if canon in seen_smiles:
            return
        seen_smiles.add(canon)
        candidates.append({
            "name"  : name,
            "smiles": canon,
            "mol"   : mol,
            "source": source,
            "linker": linker_name,
        })

    # ── ① 最重要残基の単体フラグメント ──
    if frags:
        _, name, clean = frags[0]
        _try_add(f"fragment_{name}", clean, "fragment")

    # ── リンカー決定ロジック ──
    def _linkers_for_pair(num_a: int, name_a: str,
                          num_b: int, name_b: str) -> list[dict]:
        """残基ペアに使う FEgrow リンカーエントリ群を返す (list[dict])"""
        if cbeta_coords and num_a in cbeta_coords and num_b in cbeta_coords:
            dist = float(np.linalg.norm(
                cbeta_coords[num_a] - cbeta_coords[num_b]
            ))
            ra = FRAG_REACH.get(name_a, 2.5)
            rb = FRAG_REACH.get(name_b, 2.5)
            needed = max(dist - ra - rb, 0.0)
            if needed <= 2.0:
                lks = _SHORT_LINKER_ENTRIES  # 直接連結または短鎖
            else:
                lks = _ll.select_linkers_by_dist(
                    _get_linker_db(), needed, tolerance=1.5, max_n=5,
                )
            print(f"    Cβ–Cβ({name_a}{num_a}–{name_b}{num_b}): "
                  f"{dist:.1f} Å → リンカー {[l['name'] for l in lks]}")
            return lks
        else:
            # 座標なし: 短/中/長 の3パターン
            return _FALLBACK_LINKER_ENTRIES

    # ── ② 直接連結 (距離適応リンカー、FEgrow assemble_with_linker) ──
    if n_frags >= 2:
        num_a, name_a, clean_a = frags[0]
        num_b, name_b, clean_b = frags[1]
        for lk_entry in _linkers_for_pair(num_a, name_a, num_b, name_b):
            combined = _ll.assemble_with_linker(
                clean_a, lk_entry["smiles_r1r2"], clean_b,
            )
            if combined is None:
                continue
            _try_add(
                f"direct_link_{name_a}_{name_b}_{lk_entry['name']}",
                combined, "direct_link", lk_entry["name"],
            )
        # 3残基目があれば (残基2–3) も連結 (コアSMILESによる連結)
        if n_frags >= 3:
            num_c, name_c, clean_c = frags[2]
            lks_ab = _linkers_for_pair(num_a, name_a, num_b, name_b)[:1]
            lks_bc = _linkers_for_pair(num_b, name_b, num_c, name_c)[:1]
            if lks_ab and lks_bc:
                lk_ab_core = _get_linker_core(lks_ab[0]["smiles_r1r2"]) or ""
                lk_bc_core = _get_linker_core(lks_bc[0]["smiles_r1r2"]) or ""
                combined3 = clean_a + lk_ab_core + clean_b + lk_bc_core + clean_c
                _try_add(
                    f"direct_link_{name_a}_{name_b}_{name_c}",
                    combined3, "direct_link", f"{lk_ab_core}/{lk_bc_core}",
                )

    # ── ③ レデュースドペプチド (アミド骨格 + 距離適応リンカー) ──
    if n_frags >= 2:
        # 残基間のリンカー: 短/中/長 から1つずつ (コアSMILES使用)
        for lk_idx, lk in enumerate(_FALLBACK_LINKER_CORES):
            if cbeta_coords:
                # 座標あり: 計算した最適リンカーエントリの1番目からコアを取得
                num_a, name_a, _ = frags[0]
                num_b, name_b, _ = frags[1]
                lk_entries = _linkers_for_pair(num_a, name_a, num_b, name_b)
                lk = _get_linker_core(lk_entries[0]["smiles_r1r2"]) or "" if lk_entries else ""

            backbone_parts = []
            for i, (_, name, smi) in enumerate(frags):
                side = smi + "C"   # [*] → C で Cα接続
                if i == 0:
                    backbone_parts.append(f"N{side}")
                else:
                    backbone_parts.append(f"{lk}C(=O)N{side}")
            reduced = "".join(backbone_parts)
            suffix = f"C{len(lk)}" if lk else "direct"
            _try_add(
                f"reduced_peptide_{'_'.join(n for _, n, _ in frags)}_{suffix}",
                reduced, "reduced_peptide", lk,
            )
            if cbeta_coords:
                break  # 座標あり時は1パターンで十分

    # ── ④ スキャフォールドテンプレート ──
    templates = SCAFFOLD_TEMPLATES.get(min(n_frags, 4), SCAFFOLD_TEMPLATES.get(2, []))
    for scaffold_name, scaffold_smi in templates[:2]:   # 最大2テンプレート
        if n_frags >= 2:
            smi_a = frags[0][2]
            smi_b = frags[1][2]
            smi_c = frags[2][2] if n_frags >= 3 else "C"
            composed = (scaffold_smi
                        .replace("[*:1]", smi_a)
                        .replace("[*:2]", smi_b)
                        .replace("[*:3]", smi_c))
            _try_add(
                f"{scaffold_name}_{'_'.join(n for _, n, _ in frags)}",
                composed, scaffold_name,
            )

    return candidates


def calculate_drug_likeness(mol) -> dict:
    """Lipinski Ro5 + Veber則のチェック (utils.drug_likeness に委譲)"""
    from utils.drug_likeness import calculate_drug_likeness as _calc
    return _calc(mol)


def generate_3d_conformer(mol, num_confs: int = 10) -> Optional[object]:
    """3D配座を生成して最安定構造を返す"""
    mol_h = Chem.AddHs(mol)
    params = AllChem.EmbedParameters()
    params.randomSeed = 42
    AllChem.EmbedMultipleConfs(mol_h, numConfs=num_confs, params=params)
    if mol_h.GetNumConformers() == 0:
        return None
    AllChem.MMFFOptimizeMoleculeConfs(mol_h)
    # 最低エネルギー配座を選択
    ff_props = AllChem.MMFFGetMoleculeProperties(mol_h)
    if ff_props is None:
        return mol_h
    energies = []
    for cid in range(mol_h.GetNumConformers()):
        ff = AllChem.MMFFGetMoleculeForceField(mol_h, ff_props, confId=cid)
        if ff:
            energies.append((ff.CalcEnergy(), cid))
    if energies:
        best_cid = min(energies)[1]
        # 最良配座のみ保持
        best_mol = Chem.RWMol(mol_h)
        conf = mol_h.GetConformer(best_cid)
        new_mol = Chem.RWMol(Chem.RemoveHs(mol_h))
        return new_mol
    return mol_h


def save_candidates_sdf(candidates: list[dict], output_path: str):
    """全候補をSDFに保存"""
    writer = Chem.SDWriter(output_path)
    for c in candidates:
        mol = c["mol"]
        mol_3d = generate_3d_conformer(mol)
        target = mol_3d if mol_3d else mol
        target.SetProp("_Name", c["name"])
        target.SetProp("Source", c["source"])
        target.SetProp("SMILES", c["smiles"])
        props = calculate_drug_likeness(mol)
        for k, v in props.items():
            target.SetProp(str(k), str(v))
        # 合成容易性スコアの追加
        if "SA_Score" in c:
            try:
                from synthesizability import get_synth_props_for_sdf
                for k, v in get_synth_props_for_sdf(c).items():
                    target.SetProp(k, v)
            except ImportError:
                pass
        writer.write(target)
    writer.close()
    print(f"  候補分子 SDF 保存: {output_path} ({len(candidates)}分子)")


def print_candidates(candidates: list[dict]):
    print("\n  設計された低分子候補:")
    has_sa = any("SA_Score" in c for c in candidates)
    if has_sa:
        print(f"  {'名前':<40} {'SMILES':<40} MW    LogP  Ro5 SA   QED    PAINS BRENK")
        print("  " + "-" * 130)
    else:
        print(f"  {'名前':<40} {'SMILES':<50} MW    LogP  Ro5")
        print("  " + "-" * 110)
    for c in candidates:
        props = calculate_drug_likeness(c["mol"])
        flag  = "✓" if props["Ro5"] else "✗"
        if has_sa:
            import math
            sa = c.get("SA_Score", float("nan"))
            qed = c.get("QED", float("nan"))
            sa_s  = f"{sa:.1f}" if not math.isnan(sa) else "N/A"
            qed_s = f"{qed:.3f}" if not math.isnan(qed) else " N/A "
            pains = "OK" if c.get("PAINS_OK", True) else "NG"
            brenk = "OK" if c.get("BRENK_OK", True) else "NG"
            print(f"  {c['name']:<40} {c['smiles'][:38]:<40} "
                  f"{props['MW']:<6.0f}{props['LogP']:<6.2f}{flag:<4}"
                  f"{sa_s:<5}{qed_s:<7}{pains:<6}{brenk}")
        else:
            print(f"  {c['name']:<40} {c['smiles'][:48]:<50} "
                  f"{props['MW']:<6.0f}{props['LogP']:<6.2f}{flag}")


def get_cbeta_coords_for_residues(pdb_path: str, chain_id: str,
                                   residues: list[tuple]) -> dict:
    """選択残基のCβ座標を {resnum: np.ndarray} で返す"""
    from Bio.PDB import PDBParser
    import warnings
    warnings.filterwarnings("ignore")
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("prot", pdb_path)
    coords = {}
    for model in structure:
        for chain in model:
            if chain.id != chain_id:
                continue
            for residue in chain:
                resnum = residue.id[1]
                # 選択残基のみ対象
                if not any(r[0] == resnum for r in residues):
                    continue
                atom_name = "CB" if "CB" in residue else "CA"
                if atom_name in residue:
                    coords[resnum] = np.array(residue[atom_name].get_coord())
    return coords


if __name__ == "__main__":
    import sys
    from analyze_interactions import load_structure, find_contacts, score_peptide_residues

    pdb = sys.argv[1] if len(sys.argv) > 1 else "Protein_Peptide.pdb"
    model    = load_structure(pdb)
    contacts = find_contacts(model, "A", "B")
    scores   = score_peptide_residues(contacts)
    selected = select_key_residues(scores, top_n=3)
    print(f"\n  選択残基: {[(n, name, f'{sc:.1f}') for n, name, sc in selected]}")

    # Cβ座標を取得してリンカー長を距離適応に
    cbeta = get_cbeta_coords_for_residues(pdb, "B", selected)
    print(f"  Cβ座標取得: {list(cbeta.keys())}")

    candidates = build_peptidomimetic(selected, cbeta_coords=cbeta)
    print_candidates(candidates)
    save_candidates_sdf(candidates, "candidate_ligands.sdf")
