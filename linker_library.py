"""
linker_library.py
=================
FEgrow リンカーライブラリを読み込み、距離に応じたリンカーを選択する共有モジュール。

データソース:
  - Linkers_from_FEgrow/library.sdf  : 1826分子、3D座標付き、[*:1]/[*:2]ダミー原子
  - Linkers_from_FEgrow/smiles.txt   : 2972行、タブ区切り
      col1: SMILES ([R1]/[R2]形式)
      col2: n_bonds (R1–R2間の結合数)
      col3: n_atoms (重原子数, R1/R2除く)
      col4: hba (HBAカウント)
      col5: score_r1
      col6: score_r2

キャッシュ: linker_db.json (SDF 3D距離計算済み)

主要関数:
  load_linker_db(force_rebuild)               → list[dict]
  select_linkers_by_dist(db, target, ...)     → list[dict]
  assemble_with_linker(frag_a, linker, frag_b)→ canonical SMILES | None
  build_bridge_smiles_fegrow(anc_a, lnk, anc_b) → (smiles, pharm_a_idx, pharm_b_idx) | None
  estimate_anchor_reach(smiles)               → float (Å)
"""

from pathlib import Path
import json
import numpy as np
from rdkit import Chem

BASE_DIR   = Path(__file__).parent
FEGROW_DIR = BASE_DIR / "Linkers_from_FEgrow"
SDF_PATH   = FEGROW_DIR / "library.sdf"
SMILES_TXT = FEGROW_DIR / "smiles.txt"
CACHE_PATH = BASE_DIR / "linker_db.json"

# ──────────────────────────────────────────
# 線形リンカー補完リスト (短距離 〜 ロングレンジ対応)
# FEgrow に線形エントリが少ない距離帯を補完し、一貫したカバレッジを確保する
# [R1]/[R2] 形式 SMILES (先頭 [R1] 末尾 [R2] = 線形コア抽出可能)
# ──────────────────────────────────────────
_LINEAR_SUPPLEMENT: list[dict] = [
    {"name": "none",      "smiles_r1r2": "[R1][R2]",            "dist": 0.0,  "n_bonds": 0, "n_atoms": 0, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C1",        "smiles_r1r2": "[R1]C[R2]",           "dist": 1.5,  "n_bonds": 1, "n_atoms": 1, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C2",        "smiles_r1r2": "[R1]CC[R2]",          "dist": 2.5,  "n_bonds": 2, "n_atoms": 2, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C3",        "smiles_r1r2": "[R1]CCC[R2]",         "dist": 3.8,  "n_bonds": 3, "n_atoms": 3, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C4",        "smiles_r1r2": "[R1]CCCC[R2]",        "dist": 5.0,  "n_bonds": 4, "n_atoms": 4, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C5",        "smiles_r1r2": "[R1]CCCCC[R2]",       "dist": 6.3,  "n_bonds": 5, "n_atoms": 5, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C6",        "smiles_r1r2": "[R1]CCCCCC[R2]",      "dist": 7.6,  "n_bonds": 6, "n_atoms": 6, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C7",        "smiles_r1r2": "[R1]CCCCCCC[R2]",     "dist": 9.0,  "n_bonds": 7, "n_atoms": 7, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "C8",        "smiles_r1r2": "[R1]CCCCCCCC[R2]",    "dist": 10.2, "n_bonds": 8, "n_atoms": 8, "hba": 0.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "ether_C2",  "smiles_r1r2": "[R1]COC[R2]",         "dist": 3.8,  "n_bonds": 3, "n_atoms": 3, "hba": 1.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "ether_C3",  "smiles_r1r2": "[R1]COCC[R2]",        "dist": 5.0,  "n_bonds": 4, "n_atoms": 4, "hba": 1.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "ether_C5",  "smiles_r1r2": "[R1]CCOCC[R2]",       "dist": 6.3,  "n_bonds": 5, "n_atoms": 5, "hba": 1.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "ether_C6",  "smiles_r1r2": "[R1]CCOCCCC[R2]",     "dist": 7.6,  "n_bonds": 6, "n_atoms": 6, "hba": 1.0, "score_r1": 0.0, "score_r2": 0.0},
    {"name": "amide_C3",  "smiles_r1r2": "[R1]CC(=O)N[R2]",     "dist": 4.0,  "n_bonds": 3, "n_atoms": 4, "hba": 2.0, "score_r1": 0.2, "score_r2": 0.2},
    {"name": "amide_C5",  "smiles_r1r2": "[R1]CC(=O)NCC[R2]",   "dist": 7.0,  "n_bonds": 5, "n_atoms": 6, "hba": 2.0, "score_r1": 0.2, "score_r2": 0.2},
    {"name": "amide_C7",  "smiles_r1r2": "[R1]CCC(=O)NCCC[R2]", "dist": 9.0,  "n_bonds": 7, "n_atoms": 8, "hba": 2.0, "score_r1": 0.2, "score_r2": 0.2},
]


# ══════════════════════════════════════════
# DB 構築・ロード
# ══════════════════════════════════════════

def _compute_dist_from_mol(mol) -> float | None:
    """RDKit Mol から [*:1]/[*:2] ダミー原子間の3D距離を計算する"""
    if mol is None or mol.GetNumConformers() == 0:
        return None
    conf = mol.GetConformer(0)
    dummy_idx = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 0]
    if len(dummy_idx) < 2:
        return None
    pi = np.array(conf.GetAtomPosition(dummy_idx[0]))
    pj = np.array(conf.GetAtomPosition(dummy_idx[1]))
    return round(float(np.linalg.norm(pi - pj)), 4)


def _build_db() -> list[dict]:
    """SDF と smiles.txt からリンカーDBを構築する"""
    # 1. SDF から SmileIndex → 3D距離 のマッピング
    sdf_dist: dict[int, float] = {}
    if SDF_PATH.exists():
        suppl = Chem.SDMolSupplier(str(SDF_PATH), removeHs=False)
        for mol in suppl:
            if mol is None:
                continue
            props = mol.GetPropsAsDict()
            idx_prop = props.get("SmileIndex", None)
            if idx_prop is None:
                continue
            dist = _compute_dist_from_mol(mol)
            if dist is not None:
                sdf_dist[int(idx_prop)] = dist
        print(f"  SDF 距離エントリ数: {len(sdf_dist)}")
    else:
        print(f"  警告: {SDF_PATH} が見つかりません")

    # 2. smiles.txt を読み込み DB 構築 (0-indexed 行番号 = SmileIndex)
    db: list[dict] = []
    if not SMILES_TXT.exists():
        print(f"  警告: {SMILES_TXT} が見つかりません")
        return db

    with open(SMILES_TXT, "r") as f:
        for line_num, line in enumerate(f):
            line = line.strip()
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            smiles_r1r2 = parts[0].strip()
            try:
                n_bonds  = int(parts[1])
                n_atoms  = int(parts[2])
                hba      = float(parts[3])
                score_r1 = float(parts[4]) if len(parts) > 4 else 0.0
                score_r2 = float(parts[5]) if len(parts) > 5 else 0.0
            except (ValueError, IndexError):
                continue
            db.append({
                "name":        f"fegrow_{line_num}",
                "smiles_r1r2": smiles_r1r2,
                "dist":        sdf_dist.get(line_num, None),  # None = 3D情報なし
                "n_bonds":     n_bonds,
                "n_atoms":     n_atoms,
                "hba":         hba,
                "score_r1":    score_r1,
                "score_r2":    score_r2,
            })
    return db


def load_linker_db(force_rebuild: bool = False) -> list[dict]:
    """
    FEgrow リンカー DB をロード。キャッシュ (linker_db.json) があれば使用。

    Returns:
        list[dict] with keys:
          name, smiles_r1r2, dist(Å or None), n_bonds, n_atoms, hba, score_r1, score_r2
    """
    if not force_rebuild and CACHE_PATH.exists():
        with open(CACHE_PATH, "r") as f:
            return json.load(f)

    print("FEgrow リンカーライブラリを構築中 (初回のみ)...")
    db = _build_db()
    with open(CACHE_PATH, "w") as f:
        json.dump(db, f, indent=2)
    print(f"  {len(db)} エントリをキャッシュ保存 → {CACHE_PATH}")
    return db


# ══════════════════════════════════════════
# リンカー選択
# ══════════════════════════════════════════

def select_linkers_by_dist(
    db: list[dict],
    target_dist: float,
    tolerance: float = 1.5,
    max_n: int = 5,
    max_hba: float = 4.0,
) -> list[dict]:
    """
    距離に基づいて FEgrow リンカーを選択する。

    Args:
        db         : load_linker_db() の返り値
        target_dist: 目標リンカー端-端距離 (Å)
                     = Cβ-Cβ距離 - フラグメントAリーチ - フラグメントBリーチ
        tolerance  : ±許容距離 (Å)
        max_n      : 返す最大数
        max_hba    : HBA 上限フィルタ

    Returns:
        score_r1+score_r2 降順でソートしたリスト (最大 max_n 件)
        dist が None のエントリは除外 (3D情報なし)
    """
    # ① FEgrow ライブラリから距離マッチするエントリをスコア降順で max_n 件収集
    fegrow_candidates: list[dict] = [
        e for e in db
        if e["dist"] is not None
        and e["hba"] <= max_hba
        and abs(e["dist"] - target_dist) <= tolerance
    ]
    fegrow_candidates.sort(key=lambda x: -(x["score_r1"] + x["score_r2"]))
    fegrow_result = fegrow_candidates[:max_n]

    # ② 線形補完リストから距離マッチするエントリを常に追加 (重複名除く)
    #    FEgrow が max_n を埋めていても補完は追加する (非線形FEgrowを使い切らないため)
    fegrow_names = {e["name"] for e in fegrow_result}
    supp_result: list[dict] = [
        e for e in _LINEAR_SUPPLEMENT
        if e["hba"] <= max_hba
        and abs(e["dist"] - target_dist) <= tolerance
        and e["name"] not in fegrow_names
    ]

    result = fegrow_result + supp_result

    # ③ それでも空: 全エントリ+補完から最近傍を1つ返す
    if not result:
        all_entries = [e for e in db if e["dist"] is not None] + _LINEAR_SUPPLEMENT
        if all_entries:
            filtered = [e for e in all_entries if e["hba"] <= max_hba]
            search_pool = filtered if filtered else all_entries
            best = min(search_pool, key=lambda x: abs(x["dist"] - target_dist))
            result = [best]

    return result


def estimate_anchor_reach(anchor_smiles: str) -> float:
    """
    アンカー SMILES の結合点からファーマコフォア原子までの大まかなリーチ推定 (Å)。
    簡略推定: 重原子数 × 0.9 Å
    """
    mol = Chem.MolFromSmiles(anchor_smiles)
    if mol is None:
        return 1.5
    return mol.GetNumHeavyAtoms() * 0.9


# ══════════════════════════════════════════
# SMILES 組み立て
# ══════════════════════════════════════════

def assemble_with_linker(
    frag_a_smi: str,
    linker_smi_r1r2: str,
    frag_b_smi: str,
) -> str | None:
    """
    フラグメントA + FEgrow リンカー + フラグメントB を連結して canonical SMILES を返す。

    方法:
      linker 内の [R1] → frag_a_smi、[R2] → frag_b_smi を文字列置換
      → RDKit で検証 → canonical SMILES を返す

    Args:
        frag_a_smi     : フラグメント A の SMILES
        linker_smi_r1r2: FEgrow リンカー SMILES ([R1]/[R2] 形式)
        frag_b_smi     : フラグメント B の SMILES

    Returns:
        組み立て後の canonical SMILES、失敗時は None
    """
    assembled = linker_smi_r1r2.replace("[R1]", frag_a_smi).replace("[R2]", frag_b_smi)
    mol = Chem.MolFromSmiles(assembled)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol)


def get_linker_core(smiles_r1r2: str) -> str | None:
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
    # 残りに [R1]/[R2] が含まれていれば非線形
    if "[R1]" in s or "[R2]" in s:
        return None
    return s  # "" は直接連結 (valid)


def build_bridge_smiles_linear(
    anchor_a: dict,
    linker_entry: dict,
    anchor_b: dict,
) -> tuple | None:
    """
    線形コアリンカーを使ってファーマコフォアブリッジ SMILES を構築する (文字列連結方式)。

    アンカー官能基 (carboxylic acid, phenyl 等) のように異種原子で始まる場合に
    原子価エラーを回避するため、FEgrow の [R1]/[R2] 置換は使わず
    anchor_a + core + anchor_b の単純連結でインデックスを計算する。

    非線形リンカー (core 抽出不可) は None を返す。

    Returns:
        (smiles, pharm_a_idx, pharm_b_idx) または None
    """
    core = get_linker_core(linker_entry["smiles_r1r2"])
    if core is None:
        return None  # 非線形リンカーはスキップ

    a_smi = anchor_a["smiles"]
    b_smi = anchor_b["smiles"]
    full_smiles = a_smi + core + b_smi

    mol = Chem.MolFromSmiles(full_smiles)
    if mol is None:
        return None

    a_mol = Chem.MolFromSmiles(a_smi)
    if a_mol is None:
        return None
    n_a = a_mol.GetNumHeavyAtoms()

    core_mol = Chem.MolFromSmiles(core) if core else None
    n_core = core_mol.GetNumHeavyAtoms() if core_mol else 0

    pharm_a_idx = anchor_a["pharm_offset"]
    pharm_b_idx = n_a + n_core + anchor_b["pharm_offset"]

    n_total = mol.GetNumHeavyAtoms()
    if pharm_a_idx >= n_total or pharm_b_idx >= n_total:
        return None

    return full_smiles, pharm_a_idx, pharm_b_idx


def build_bridge_smiles_fegrow(
    anchor_a: dict,
    linker_entry: dict,
    anchor_b: dict,
) -> tuple | None:
    """
    FEgrow リンカーを使ってファーマコフォアブリッジ SMILES を構築し、
    ファーマコフォア原子インデックスを atom map 番号で追跡する。

    Args:
        anchor_a    : ANCHOR_A エントリ (smiles, pharm_offset, name, ...)
        linker_entry: FEgrow リンカーエントリ (smiles_r1r2, ...)
        anchor_b    : ANCHOR_B エントリ (smiles, pharm_offset, name, ...)

    Returns:
        (canonical_smiles, pharm_a_idx, pharm_b_idx) または None
    """
    a_smi = anchor_a["smiles"]
    b_smi = anchor_b["smiles"]
    linker_smi = linker_entry["smiles_r1r2"]

    a_mol = Chem.MolFromSmiles(a_smi)
    b_mol = Chem.MolFromSmiles(b_smi)
    if a_mol is None or b_mol is None:
        return None

    pharm_a_off = anchor_a["pharm_offset"]
    pharm_b_off = anchor_b["pharm_offset"]

    if pharm_a_off >= a_mol.GetNumHeavyAtoms():
        return None
    if pharm_b_off >= b_mol.GetNumHeavyAtoms():
        return None

    # ファーマコフォア原子に atom map 番号を付与
    a_rw = Chem.RWMol(a_mol)
    b_rw = Chem.RWMol(b_mol)
    a_rw.GetAtomWithIdx(pharm_a_off).SetAtomMapNum(91)  # 91/92 は衝突しにくい値
    b_rw.GetAtomWithIdx(pharm_b_off).SetAtomMapNum(92)

    a_smi_mapped = Chem.MolToSmiles(a_rw)
    b_smi_mapped = Chem.MolToSmiles(b_rw)

    # リンカーを介して組み立て
    assembled_smi = linker_smi.replace("[R1]", a_smi_mapped).replace("[R2]", b_smi_mapped)
    mol_mapped = Chem.MolFromSmiles(assembled_smi)
    if mol_mapped is None:
        return None

    # atom map 番号からインデックスを取得
    pharm_a_idx = -1
    pharm_b_idx = -1
    for atom in mol_mapped.GetAtoms():
        mn = atom.GetAtomMapNum()
        if mn == 91:
            pharm_a_idx = atom.GetIdx()
        elif mn == 92:
            pharm_b_idx = atom.GetIdx()

    if pharm_a_idx == -1 or pharm_b_idx == -1:
        return None

    # atom map 番号を除去して最終 SMILES を取得
    clean_mol = Chem.RWMol(mol_mapped)
    for atom in clean_mol.GetAtoms():
        atom.SetAtomMapNum(0)
    final_smi = Chem.MolToSmiles(clean_mol.GetMol())

    return final_smi, pharm_a_idx, pharm_b_idx


# ══════════════════════════════════════════
# CLI (テスト用)
# ══════════════════════════════════════════

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="FEgrow リンカーライブラリ管理")
    parser.add_argument("--rebuild",  action="store_true", help="キャッシュを再構築")
    parser.add_argument("--dist",     type=float, default=None, help="距離でリンカー検索")
    parser.add_argument("--tol",      type=float, default=1.5,  help="許容誤差 (Å)")
    parser.add_argument("--max-n",    type=int,   default=10,   help="最大表示数")
    parser.add_argument("--stats",    action="store_true",      help="距離分布を表示")
    args = parser.parse_args()

    db = load_linker_db(force_rebuild=args.rebuild)
    print(f"\nDB エントリ数: {len(db)}")

    dists = [e["dist"] for e in db if e["dist"] is not None]
    print(f"3D距離あり: {len(dists)} / {len(db)}")
    if dists and args.stats:
        import statistics
        print(f"  距離範囲: {min(dists):.2f} – {max(dists):.2f} Å")
        print(f"  平均: {statistics.mean(dists):.2f} Å, 中央値: {statistics.median(dists):.2f} Å")
        bins = [(i, i+1) for i in range(0, 12)]
        for lo, hi in bins:
            cnt = sum(1 for d in dists if lo <= d < hi)
            if cnt:
                print(f"  {lo}–{hi} Å: {cnt}")

    if args.dist is not None:
        results = select_linkers_by_dist(db, args.dist, args.tol, args.max_n)
        print(f"\n目標距離 {args.dist:.1f} Å ± {args.tol} Å → {len(results)} 件")
        for r in results:
            dist_str = f"{r['dist']:.2f}" if r["dist"] is not None else "N/A"
            print(f"  {r['name']:20s}  dist={dist_str:5s}  hba={r['hba']:.1f}"
                  f"  score={r['score_r1']+r['score_r2']:+.3f}  SMILES={r['smiles_r1r2']}")
