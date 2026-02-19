# Peptide → Small Molecule Conversion Platform

ペプチド-タンパク質複合体の**結合相互作用面を保持**しながら、ペプチドを低分子リガンドに変換し、sminaでドッキング検証するパイプライン。

---

## 目次

1. [概要](#概要)
2. [ファイル構成](#ファイル構成)
3. [アルゴリズム詳細](#アルゴリズム詳細)
   - [STEP 1: 距離ベース相互作用解析](#step-1-距離ベース相互作用解析)
   - [STEP 1b: ΔSASA解析](#step-1b-δsasa解析)
   - [STEP 1c: PRODIGY結合親和性予測](#step-1c-prodigy結合親和性予測)
   - [STEP 2: ファーマコフォア抽出](#step-2-ファーマコフォア抽出)
   - [STEP 3: 低分子設計](#step-3-低分子設計)
   - [STEP 4: 可視化](#step-4-可視化)
   - [STEP 5: sminaドッキング](#step-5-sminaドッキング)
   - [STEP 6b: sminaスコア項分解](#step-6b-sminaスコア項分解)
   - [スコア統合](#スコア統合)
   - [STEP 7: ファーマコフォアブリッジ分子生成](#step-7-ファーマコフォアブリッジ分子生成)
   - [STEP 8: 環状ペプチド vs 低分子 比較 (smina)](#step-8-環状ペプチド-vs-低分子-比較)
   - [STEP 9: AutoDock CrankPep 環状ペプチドドッキング](#step-9-autodock-crankpep-環状ペプチドドッキング)
   - [STEP 10: AutoDockFR 低分子ドッキング](#step-10-autodockfr-低分子ドッキング)
4. [手法の比較](#手法の比較)
5. [使い方](#使い方)
6. [出力ファイル](#出力ファイル)
7. [最終選抜結果 (Result_Best)](#最終選抜結果-result_best)
8. [依存関係](#依存関係)

---

## 概要

```
Protein_Peptide.pdb
        │
        ▼
┌─────────────────────────────────────────────────────┐
│  STEP 1   距離ベース相互作用解析 (analyze_interactions)│
│           接触検出 → H結合/疎水性/静電 分類           │
│           → 残基別スコアリング (加重和)               │
├─────────────────────────────────────────────────────┤
│  STEP 1b  ΔSASA解析 (analyze_sasa)                  │
│           複合体 vs 単体 の表面積差分                 │
│           → 界面埋没面積 (BSA) / コア残基同定         │
├─────────────────────────────────────────────────────┤
│  STEP 1c  PRODIGY結合親和性予測 (analyze_prodigy)    │
│           Cβ/Cα接触タイプ(CC/CP/CN/PP/PN/NN)カウント │
│           → 線形回帰で ΔG・Kd 予測                   │
│           → 距離スコア × ΔSASA の統合スコア生成       │
└───────────────────────┬─────────────────────────────┘
                        │ 統合残基ランキング
              ┌─────────┴──────────┐
              ▼                    ▼
┌──────────────────────┐  ┌──────────────────────────────────┐
│  STEP 2              │  │  STEP 7   ファーマコフォアブリッジ │
│                      │  │                                  │
│  ファーマコフォア抽出 │  │  (pharmacophore_bridge)          │
│  HBA/HBD/HYD/ARO/    │  │  タンパク質側ポケット残基2点を指定│
│  POS/NEG             │  │  → Cβ–Cβ距離を目標に            │
│  → 3D座標特徴点      │  │  → アンカー + リンカー + アンカー│
└──────────┬───────────┘  │  → RDKit距離幾何学 (DG法) で     │
           │              │    3D配座を距離拘束付き生成       │
           ▼              │  → MMFF最適化 → smina ドッキング │
┌──────────────────────┐  └──────────────────────────────────┘
│  STEP 3              │
│  低分子設計          │
│  (ペプチドミメティクス)│
│  側鎖フラグメント抽出 │
│  → 骨格連結 (4戦略)  │
│  → 3D配座生成 → SDF  │
└──────────┬───────────┘
           │ candidate_ligands.sdf
           ▼
┌─────────────────────────────────────────────────────┐
│  STEP 6  smina ドッキング (dock_with_smina)          │
│          autobox_ligand でボックス自動定義            │
│          → アフィニティスコア順ランキング             │
├─────────────────────────────────────────────────────┤
│  STEP 6b smina --score_only リスコアリング           │
│          gauss / hydrophobic / Hbond / repulsion     │
│          → スコア項分解グラフ (MM-GBSA近似)           │
└─────────────────┬───────────────────────────────────┘
                  │ 全結果 (CSV) を集約
                  ▼
┌─────────────────────────────────────────────────────┐
│  STEP 8  環状ペプチド vs 低分子 比較 (smina)          │
│  (compare_cyclic_peptide)                           │
│  cyclo-(配列) SMILES 自動生成                        │
│  → ETKDG + MMFF で 3D 配座生成                      │
│  → smina --score_only (線形ペプチド, ネイティブポーズ)│
│  → smina docking (環状ペプチド de novo)              │
│  → 全設計低分子と比較グラフ + PRODIGY ΔG 参照線      │
├─────────────────────────────────────────────────────┤
│  STEP 9  AutoDock CrankPep 環状ペプチドドッキング    │
│  (dock_cyclic_adcp)                                 │
│  receptor.pdb → .trg (agfr)                        │
│  → adcp -cyc (MC + implicit 溶媒)                  │
│  → 環状ペプチド Affinity (kcal/mol)                 │
│  → smina 低分子結果と並列比較 + LE グラフ            │
├─────────────────────────────────────────────────────┤
│  STEP 10 AutoDockFR 低分子ドッキング                 │
│  (dock_smol_adfr)                                   │
│  Result_Best 上位分子 SDF → PDBQT (obabel)         │
│  → adfr (GA / AutoDock 力場 / 静電+脱溶媒和)        │
│  → .trg 共用 (ADCP と同一受容体ファイル)             │
│  → ADFR vs smina / ADFR vs ADCP / LE 3 メソッド比較│
└─────────────────────────────────────────────────────┘
```

入力PDBの例（本プロジェクト）:
- **Chain A** : タンパク質受容体（115残基）
- **Chain B** : ペプチドリガンド（10残基 `GEVDGWATPD`）

---

## ファイル構成

```
Peptide_to_small_molecules/
├── pipeline.py                 # メインパイプライン（全STEP統合）
├── analyze_interactions.py     # STEP 1:  距離ベース相互作用解析
├── analyze_sasa.py             # STEP 1b: ΔSASA解析
├── analyze_prodigy.py          # STEP 1c: PRODIGY結合親和性予測
├── extract_pharmacophore.py    # STEP 2:  ファーマコフォア抽出
├── design_small_molecule.py    # STEP 3:  低分子設計（ペプチドミメティクス）
├── generate_pocket_molecules.py# STEP 3b: ポケット相補分子生成（BRICS）
├── visualize.py                # STEP 4:  可視化
├── dock_with_smina.py          # STEP 6/6b: ドッキング + rescoring
├── pharmacophore_bridge.py     # STEP 7:  ファーマコフォアブリッジ分子生成
├── compare_cyclic_peptide.py   # STEP 8:  環状ペプチドと設計低分子の比較 (smina)
├── dock_cyclic_adcp.py         # STEP 9:  AutoDock CrankPep 環状ペプチドドッキング + smina 比較
├── dock_smol_adfr.py           # STEP 10: AutoDockFR 低分子ドッキング + ADCP/smina 比較
├── adcpsuite_micromamba.sh     # ADCP/ADFR 環境セットアップスクリプト (micromamba)
├── collect_best.py             # 最終候補選抜・SDF集約スクリプト → Result_Best/
├── linker_library.py           # 共有モジュール: FEgrowリンカーライブラリ管理
├── plot_utils.py               # matplotlib 日本語フォント設定
├── smina.osx.12                # sminaバイナリ
├── Protein_Peptide.pdb         # 入力構造
├── linker_db.json              # FEgrowリンカーDB（初回実行時に自動生成・キャッシュ）
├── Linkers_from_FEgrow/        # FEgrowリンカーデータ
│   ├── library.sdf             # 1826分子、[*:1]/[*:2]ダミー原子付き3D構造
│   └── smiles.txt              # 2972行、[R1]/[R2]形式SMILES + 各種特性
├── Result_Best/                # 最終選抜候補 (collect_best.py で自動生成)
│   ├── 01_score*_<name>.sdf   # 選抜分子ドッキングポーズ (PyMOL/VMD で可視化可)
│   ├── ...                    # 02〜N 番まで同様
│   ├── summary.csv            # 全選抜分子の詳細データ
│   └── README.md              # 選抜基準・分子一覧
└── results/                    # 出力ディレクトリ
    ├── candidate_ligands.sdf
    ├── pharmacophore.csv
    ├── pharmacophore.pml       # PyMOL可視化スクリプト
    ├── contacts.json
    ├── residue_scores.png
    ├── interaction_map.png
    ├── pharmacophore_3d.png
    ├── docking/
    │   ├── receptor.pdb
    │   ├── peptide_ref.pdb
    │   ├── docked_*.sdf
    │   ├── docking_results.csv
    │   └── docking_scores.png
    ├── bridge/                 # STEP 7 出力
    │   ├── bridge_*_results.csv
    │   ├── bridge_*_analysis.png
    │   ├── geometry_*.png
    │   ├── input_bridge_*.sdf
    │   └── docked_bridge_*.sdf
    └── adcp_docking/           # STEP 9 出力 (AutoDock CrankPep)
        ├── receptor.pdb
        ├── receptor_rec.pdbqt
        ├── peptide_ref.pdbqt
        ├── receptor.trg
        ├── gevdgwatpd/
        │   ├── result_gevdgwatpd_summary.dlg
        │   └── result_gevdgwatpd_out.pdb
        ├── adcp_vs_smina_comparison.png
        ├── adcp_vs_smina_le_comparison.png
        └── adcp_comparison_report.json
```

---

## アルゴリズム詳細

### STEP 1: 距離ベース相互作用解析

**実装**: `analyze_interactions.py`

#### 1-1. 接触検出

BioPythonの `NeighborSearch` を使い、ペプチド（Chain B）の各重原子から半径 **4.5 Å** 以内にあるタンパク質（Chain A）の重原子をすべて列挙する。

```
for atom_p in peptide_heavy_atoms:
    neighbors = NeighborSearch.search(atom_p.coord, cutoff=4.5Å)
    contacts += [n for n in neighbors if n.chain == protein_chain]
```

水素原子（H）は除外する。

#### 1-2. 相互作用タイプ分類

各接触ペアを以下のルールで分類する：

| タイプ | 条件 |
|--------|------|
| **H結合 (hbond)** | N/O … N/O かつ距離 ≤ 3.5 Å |
| **疎水性 (hydrophobic)** | C … C かつ距離 ≤ 4.5 Å かつペプチド残基が疎水性AA |
| **静電気 (electrostatic)** | 両側が荷電残基 (ASP/GLU/LYS/ARG) かつ距離 ≤ 6.0 Å |
| **ファンデルワールス** | C … C かつ非疎水性残基 |
| **その他** | 上記に該当しない接触 |

疎水性アミノ酸: `ALA, VAL, ILE, LEU, MET, PHE, TRP, PRO, TYR`

#### 1-3. 残基重要度スコアリング

各ペプチド残基に対して、接触タイプ別の加重和でスコアを算出する：

```
Score(残基) = Σ weight(接触タイプ)

weight:
  hbond          = 3.0  (特異的・方向性あり → 高重要)
  electrostatic  = 2.0
  hydrophobic    = 1.0
  van_der_waals  = 0.5
  other          = 0.3
```

スコアが高い残基を「ホットスポット残基」として低分子設計に使用する。

本PDBでの結果:
```
TRP6 : 56.2点  ← 芳香族疎水性接触が多数
ALA7 : 25.9点
ASP4 : 11.8点  ← 静電相互作用
```

---

### STEP 1b: ΔSASA解析

**実装**: `analyze_sasa.py`

BioPythonの `ShrakeRupley` アルゴリズムを用いて、複合体・タンパク質単体・ペプチド単体それぞれのSASAを計算し、残基ごとの埋没量を算出する。

#### 計算フロー

```
①  SASA(複合体)     を全残基で計算
②  SASA(タンパク質単体) を計算
③  SASA(ペプチド単体)   を計算

ΔSASA(残基i) = SASA_単体(i) - SASA_複合体(i)
             ↑ 正の値 = 複合体形成で埋没 = 界面残基
```

#### 探針半径 (Probe Radius)
水分子を模した球 (r = **1.4 Å**) を分子表面に沿って転がし、接触できる面積を積算する（Lee & Richards 法）。

#### 界面残基の分類
| 閾値 | 分類 |
|------|------|
| ΔSASA ≥ 1.0 Å² | 界面残基 (interface) |
| ΔSASA ≥ 10.0 Å² | コア界面残基 (core) |

#### Buried Surface Area (BSA)
複合体全体の埋没面積 = タンパク質側ΔSASA + ペプチド側ΔSASA の合計。
```
本PDB結果: BSA = 1111.4 Å²  (ペプチド-タンパク質としては中程度)
典型値: 抗体-抗原 1500〜3000 Å², ペプチド-タンパク質 500〜1500 Å²
```

最重要残基:
```
TRP6  208 Å²   (コア界面 — インドール環が完全に埋没)
PRO9  103 Å²   (コア界面)
ALA7   99 Å²   (コア界面)
```

---

### STEP 1c: PRODIGY結合親和性予測

**実装**: `analyze_prodigy.py` — Vangone & Bonvin (eLife 2015) の再実装

#### 残基タイプ分類

| タイプ | 記号 | 残基 |
|--------|------|------|
| 荷電 | **C** | ASP, GLU, LYS, ARG, HIS |
| 極性 | **P** | SER, THR, ASN, GLN, TYR, TRP, CYS |
| 非極性 | **N** | ALA, VAL, ILE, LEU, MET, PHE, PRO, GLY |

#### 界面接触カウント (Cβ/Cα 間距離 ≤ 5.5 Å)

距離ベースの接触検出(4.5 Å, 重原子)より長いカットオフを用いるのは、Cβ/Cαは側鎖の「代表点」であり、実際の側鎖原子はより遠くまで伸びているため。

```
6タイプ: IC_CC, IC_CP, IC_CN, IC_PP, IC_PN, IC_NN
```

#### ΔG 予測式 (線形回帰)

```
ΔG = -0.09459 × IC_CC
   + -0.10007 × IC_CP
   +  0.19577 × IC_CN    ← 荷電-非極性は不利
   + -0.22671 × IC_PP
   +  0.18681 × IC_PN    ← 極性-非極性は不利
   +  0.34282 × IC_NN    ← 非極性同士も意外と不利 (脱溶媒コスト)
   +  0.01384 × NIS_P    ← 非界面極性残基 (エントロピーコスト)
   +  0.02064 × NIS_C
   + -15.9433            ← 切片
```

NIS (Non-Interface Surface) = 表面に露出しているが界面に参加していない残基数。

#### Kd 計算

```python
Kd = exp(ΔG / (R × T))   # R=1.987×10⁻³ kcal/(mol·K), T=298.15 K
```

本PDB結果: `ΔG = -14.7 kcal/mol, Kd ≈ 16 pM` (非常に強い結合を示唆)

---

### スコア統合

STEP 1（距離ベース）と STEP 1b（ΔSASA）のスコアを正規化して加重平均で統合し、より頑健な残基重要度ランキングを生成する。

```
combined(残基i) = 0.6 × (dist_score / max_dist_score)
               + 0.4 × (ΔSASA / max_ΔSASA)
```

統合スコアを STEP 3 の低分子設計に使用することで、単純な接触数だけでなく実際の「埋没量」も考慮した設計が可能になる。

---

### STEP 2: ファーマコフォア抽出

**実装**: `extract_pharmacophore.py`

STEP 1でスコア > 0 の残基について、側鎖原子の3D座標から薬理活性特徴点を抽出する。

#### 特徴タイプと対応原子

| タイプ | 記号 | 対象原子名の例 |
|--------|------|----------------|
| 水素結合受容体 | **HBA** | O, OD1, OD2, OE1, OE2, OG, OH |
| 水素結合供与体 | **HBD** | N, ND1, NE1, NE2, NZ, OG, OH |
| 疎水性 | **HYD** | C系側鎖原子 |
| 芳香族 | **ARO** | 芳香環重心 |
| 正電荷 | **POS** | NZ (LYS), NH1/NH2 (ARG) |
| 負電荷 | **NEG** | OD1/OD2 (ASP), OE1/OE2 (GLU) |

複数原子が関与する場合（例: 芳香環）は**重心座標**を使用する：

```python
center = mean([atom.coord for atom in ring_atoms])
```

出力はCSV（3D座標）とPyMOL描画スクリプト（.pml）。

---

### STEP 3: 低分子設計

**実装**: `design_small_molecule.py`

ペプチドの「非ペプチド等価体」を生成するペプチドミメティクスアプローチ。

#### 3-1. 残基 → 側鎖フラグメント変換

重要残基の側鎖をRDKit SMILES フラグメントにマッピングする（`[*]` が結合点）：

```
TRP → c1ccc2[nH]ccc2c1C[*]   (インドール環)
ALA → C[*]                    (メチル基)
ASP → OC(=O)C[*]              (カルボキシル)
```

#### 3-2. 骨格連結戦略（4種）

選択した上位残基（デフォルト3残基）を以下の方法で連結して候補分子を生成する。
`cbeta_coords` を渡すと残基間Cβ–Cβ距離からリンカー長を自動選択し、渡さない場合は短(C) / 中(CCCC) / 長(CCCCCC) の3パターンを生成する。

**① フラグメント単体**
最重要残基の側鎖をそのまま低分子として使用。

```
TRP側鎖 → Cc1cccc2[nH]ccc12 (MW=131)
```

**② 距離適応リンカー連結 (Direct Link)**
Cβ–Cβ距離とフラグメントの「リーチ」（Cαから官能基端までの推定距離）から必要なリンカー長を計算し、適切な直鎖アルキルリンカーで連結する。

```
必要リンカー長 = max(Cβ–Cβ距離 − reach_A − reach_B, 0)

# 例: TRP6–ALA7 (Cβ–Cβ = 5.8 Å)
reach_TRP = 5.5 Å, reach_ALA = 1.3 Å
必要長 = max(5.8 − 5.5 − 1.3, 0) = 0 Å → 直結 "" or "C"

# 例: ALA7–ASP4 (Cβ–Cβ = 10.4 Å)
reach_ALA = 1.3 Å, reach_ASP = 3.8 Å
必要長 = 10.4 − 1.3 − 3.8 = 5.3 Å → "CCC"〜"CCCCC" (±1.5 Å 許容)
```

リンカー選択は `linker_library.py` の `select_linkers_by_dist()` で行われる。FEgrowライブラリ (2972エントリ) を距離でフィルタリングし、補完リスト `_LINEAR_SUPPLEMENT` を常に追加する。

**`_LINEAR_SUPPLEMENT`（線形リンカー補完リスト, 一部）**:

| リンカー名 | `smiles_r1r2` | 伸展長 (Å) |
|-----------|--------------|-----------|
| none      | `[R1][R2]`          | 0.0  |
| C1        | `[R1]C[R2]`         | 1.5  |
| C2        | `[R1]CC[R2]`        | 2.5  |
| C3        | `[R1]CCC[R2]`       | 3.8  |
| C4        | `[R1]CCCC[R2]`      | 5.0  |
| C5        | `[R1]CCCCC[R2]`     | 6.3  |
| C6        | `[R1]CCCCCC[R2]`    | 7.6  |
| C7        | `[R1]CCCCCCC[R2]`   | 9.0  |
| ether_C2  | `[R1]COC[R2]`       | 3.8  |
| ether_C5  | `[R1]CCOCC[R2]`     | 6.3  |
| amide_C3  | `[R1]CC(=O)N[R2]`   | 4.0  |
| amide_C5  | `[R1]CC(=O)NCC[R2]` | 7.0  |

FEgrowリンカーは `assemble_with_linker()` で `[R1]→フラグメントA`, `[R2]→フラグメントB` に置換して組み立てる。FEgrowに線形エントリが少ない距離帯では `_LINEAR_SUPPLEMENT` が自動的に補完する（最大距離10.2 Åまで対応）。

**③ レデュースドペプチド (Reduced Peptide)**
アミド結合（-CO-NH-）のバックボーンを残しながら、残基間リンカー長も距離適応で調整する。N末端にアミノ基を配置して非ペプチド的性質を持たせる。

```
NH2-TRP_side - [linker] - C(=O)NH-ALA_side - [linker] - C(=O)NH-ASP_side
```

**④ スキャフォールドテンプレート**
ベンゼン環やピペラジンなどの剛直な骨格に側鎖を配置（最大2テンプレート）。

```
ベンゼン三置換:  c1c([R1])cc([R2])cc1[R3]
ピペラジン:     C1CN(CC([R1])[R2])CC1
エチレン:       [R1]CC[R2]
```

#### 3-3. 3D配座生成

RDKitの `EmbedMultipleConfs` + `MMFFOptimizeMoleculeConfs` で10配座を生成し、MMFF94力場エネルギーが最小の配座を選択してSDF出力する。

#### 3-4. ドラッグライクネス評価

Lipinskiのルールオブファイブ (Ro5) と Veber則でフィルタリング：

```
Ro5 : MW ≤ 500, LogP ≤ 5, HBD ≤ 5, HBA ≤ 10
Veber: 回転可能結合 ≤ 10, PSA ≤ 140 Å²
```

---

### STEP 4: 可視化

**実装**: `visualize.py`

| 出力ファイル | 内容 |
|---|---|
| `residue_scores.png` | ペプチド残基ごとの重要度スコア棒グラフ（上位3残基を赤/橙/黄でハイライト） |
| `interaction_map.png` | タンパク質残基 × ペプチド残基の接触数ヒートマップ（上位20残基を表示） |
| `pharmacophore_3d.png` | ファーマコフォア特徴点の3D散布図（特徴タイプ別に色分け） |
| `docking_scores.png` | 候補分子のドッキングスコア水平棒グラフ |

---

### STEP 5: sminaドッキング

**実装**: `dock_with_smina.py`

#### 5-1. 受容体・参照リガンドの準備

PDB から Chain A（受容体）と Chain B（参照ペプチド）を別ファイルに分割する。

```python
# chain_id列 (PDB列21) でフィルタリング
if line[21] == chain_id:
    keep(line)
```

#### 5-2. ドッキングボックス自動定義

`--autobox_ligand` オプションに参照ペプチド（Chain B）を渡すことで、**ペプチドの結合部位を包含するボックスをsminaが自動計算**する。

```
smina --autobox_ligand peptide_ref.pdb --autobox_add 4.0
```

`autobox_add 4.0` でペプチド座標範囲の各辺に +4 Å のバッファを追加（最小20 Å確保）。

手動でのボックス計算も実装済み（バックアップ用）：

```python
center = (max_coord + min_coord) / 2
size   = (max_coord - min_coord) + 2 * buffer
```

#### 5-3. ドッキング実行

候補分子をSDF1分子ずつ個別にドッキングすることで、スコアを分子ごとに追跡できる。

```bash
smina \
  --receptor  receptor.pdb        \
  --ligand    input_mol.sdf       \
  --autobox_ligand peptide_ref.pdb\
  --autobox_add 4.0               \
  --out       docked_mol.sdf      \
  --exhaustiveness 8              \
  --num_modes 9                   \
  --energy_range 5
```

#### 5-4. スコア解析

sminaログからアフィニティスコアを正規表現で抽出：

```
パターン: ^\s*\d+\s+([-\d.]+)\s+[\d.]+\s+[\d.]+
                    ^^^^^^^^
                    affinity (kcal/mol)
```

各分子のベストスコアでランキングし、CSV・JSON・グラフを出力する。

---

### STEP 6b: sminaスコア項分解

**実装**: `dock_with_smina.py` の `rescore_with_smina()`

ドッキングで得られたポーズに対して `smina --score_only` を適用し、Vina スコア関数の各項を個別に取得する。これにより「なぜこのポーズが良いスコアか」を定量的に分解できる（MM-GBSAエネルギー分解の簡易版）。

#### smina のスコア関数

```
ΔG_Vina = w₁·gauss₁ + w₂·gauss₂ + w₃·repulsion
        + w₄·hydrophobic + w₅·Hbond + num_tors_div

w₁ = -0.035579  (gauss1: 短距離vdW引力)
w₂ = -0.005156  (gauss2: 長距離vdW引力)
w₃ =  0.840245  (repulsion: 立体反発)
w₄ = -0.035069  (hydrophobic: 疎水性接触)
w₅ = -0.587439  (Hbond: 水素結合)
```

smina の `--score_only` 出力フォーマット：

```
## Name  gauss(o=0,...) gauss(o=3,...) repulsion(...) hydrophobic(...) h_bond(...) num_tors_div
Affinity: -7.154 (kcal/mol)
Term values, before weighting:
## mol_name  raw₁  raw₂  raw₃  raw₄  raw₅  raw₆
```

重みつき寄与 = `重み × raw値` で各項の貢献を計算する。

本PDB結果例 (`reduced_peptide_TRP_ALA_PRO`、ベストスコア分子):
```
Affinity:    -7.15 kcal/mol
H結合:       -2.22 kcal/mol  ← 最大の貢献
疎水性:      -1.54 kcal/mol
立体反発:    +2.39 kcal/mol  (不利な要因)
```

---

### 共有モジュール: linker_library.py

`linker_library.py` は `design_small_molecule.py`（STEP 3）と `pharmacophore_bridge.py`（STEP 7）の両方から使用される共有リンカー管理モジュール。

| 関数 | 説明 |
|------|------|
| `load_linker_db(force_rebuild)` | SDF+smiles.txt を読み込み `linker_db.json` にキャッシュ |
| `select_linkers_by_dist(db, target, ...)` | 距離でフィルタしたFEgrowエントリ + `_LINEAR_SUPPLEMENT` を返す |
| `get_linker_core(smiles_r1r2)` | `[R1]xxx[R2]` → `xxx`（非線形は `None`）|
| `assemble_with_linker(frag_a, lnk, frag_b)` | `[R1]→frag_a`, `[R2]→frag_b` の文字列置換（STEP 3用）|
| `build_bridge_smiles_linear(anchor_a, lnk, anchor_b)` | 線形コア連結 + ファーマコフォア原子索引計算（STEP 7用）|
| `estimate_anchor_reach(smiles)` | アンカーの伸展距離推定: `n_重原子 × 0.9 Å` |

`linker_db.json` は初回実行時に自動生成される（`Linkers_from_FEgrow/` フォルダが存在する場合）。

---

### STEP 7: ファーマコフォアブリッジ分子生成

**実装**: `pharmacophore_bridge.py`（共有モジュール `linker_library.py` を使用）

ポケット内の2つのタンパク質残基を指定し、そのCβ–Cβ距離を「橋渡し距離」の目標として、**両残基と同時に相互作用できる低分子**を3D構造拘束付きで設計するモジュール。

**従来のペプチドミメティクス（STEP 3）との差異**:

| 手法 | 設計基準 | 3D拘束 |
|------|----------|--------|
| STEP 3 (ペプチドミメティクス) | ペプチド側鎖を模倣 | なし (ランダム配座生成) |
| STEP 7 (ファーマコフォアブリッジ) | タンパク質ポケット形状に合わせて設計 | あり (Cβ–Cβ距離制約) |

---

#### 7-1. ポケット残基ペアと目標距離の算出

PDBからタンパク質側の2残基（ユーザー指定）のCβ座標を読み取り、Cβ–Cβ距離を目標距離 $d_{\text{target}}$ とする。

```
GLYの場合はCβが存在しないためCαを使用する。
```

本PDBで解析済みのポケットペア（タンパク質側, Chain A）:

| ペア | Cβ–Cβ距離 | タイプ |
|------|-----------|--------|
| LYS48 — LEU50 | 7.3 Å | POS–HYD |
| LYS46 — TYR49 | 6.6 Å | POS–ARO |
| ILE21 — TYR25 | 6.1 Å | HYD–ARO |
| ARG28 — TYR49 | 9.1 Å | POS–ARO |
| LYS46 — LEU50 | 10.0 Å | POS–HYD |

---

#### 7-2. 相補的アンカー基の選択

各残基タイプに対して「タンパク質残基と相互作用できる」相補的な官能基（アンカー）を選択する。選択ルールはSTEP 3bの相補性ロジックを踏襲。

| タンパク質残基タイプ | アンカー基候補 | 相互作用様式 |
|---------------------|--------------|-------------|
| POS (LYS/ARG) | カルボン酸, スルホンアミド, テトラゾール | 塩橋 / H結合 |
| NEG (ASP/GLU) | 一級アミン, グアニジン | 塩橋 |
| ARO (PHE/TYR/TRP) | フェニル, ナフタレン, インドール | π–πスタッキング |
| HYD (LEU/ILE/VAL) | シクロヘキシル, フェニル, イソプロピル | 疎水性接触 |
| POL (SER/THR/ASN) | ヒドロキシル, アミド | 水素結合 |

---

#### 7-3. リンカー長の自動選択

アンカーA・アンカーBの「リーチ（伸展距離）」をそれぞれ個別に計算し、必要なリンカー長（`linker_target`）を算出する。この計算はアンカーペアごとに行われる。

```python
# アンカーごとのリーチ推定: n_重原子 × 0.9 Å
a_reach = estimate_anchor_reach(anchor_a["smiles"])
b_reach = estimate_anchor_reach(anchor_b["smiles"])

# 必要リンカー長 = Cβ–Cβ目標距離 - 両アンカーのリーチ合計
linker_target = max(0.0, target_dist - a_reach - b_reach)

# リンカー選択 (tolerance = ±2.5 Å)
linkers = select_linkers(linker_target, tolerance=2.5)
```

**`select_linkers_by_dist()` のアルゴリズム**:

```
① FEgrow ライブラリ (Linkers_from_FEgrow/library.sdf + smiles.txt)
   距離 = SDF中の [*:1]–[*:2] ダミー原子間3D距離
   linker_target ± tolerance の範囲でフィルタ
   スコア降順で最大 max_n=5 件を選択

② _LINEAR_SUPPLEMENT から距離マッチ分を常に追加
   (FEgrow の件数に関わらず追加 — max_n に縛られない)

③ 候補が空の場合: 全エントリの中から最近傍1件を返す
```

FEgrow ライブラリのスペック:

| 項目 | 値 |
|------|---|
| 総エントリ数 | 2972 (smiles.txt) |
| 3D距離付き | 1825 (library.sdf) |
| 距離範囲 | 2.32 – 6.47 Å (平均 3.93 Å) |
| 線形コア抽出可能 | 約47% (864件) |
| `_LINEAR_SUPPLEMENT` でカバー | 0.0 – 10.2 Å |

---

#### 7-4. SMILES組み立てとファーマコフォア原子索引の追跡

`build_bridge_smiles_linear()` (`linker_library.py`) を使い、リンカーエントリから線形コアを抽出して文字列結合でSMILESを組み立てる：

```python
# リンカーエントリ {"smiles_r1r2": "[R1]CCCC[R2]", ...} から線形コアを取得
core = get_linker_core(linker_entry["smiles_r1r2"])  # "[R1]CCCC[R2]" → "CCCC"
# コアが取得できない (非線形FEgrowリンカー等) → None を返してスキップ
if core is None:
    return None

full_SMILES = AnchorA_SMILES + core + AnchorB_SMILES
```

`get_linker_core()` は `[R1]xxx[R2]` 形式から `xxx` を取り出す。環への付加・分岐など非線形構造の場合は `None` を返し、そのリンカーエントリは自動的にスキップされる。

各パーツの重原子数 $n_A, n_{\text{core}}, n_B$ を用いて、ファーマコフォア原子の索引を決定する：

```
pharm_A_idx = AnchorA内のオフセット (残基タイプ依存)
pharm_B_idx = n_A + n_core + AnchorB内のオフセット
```

例: LYS(POS) → LEU(HYD), リンカー=C4

```
AnchorA  = "OC(=O)"        # COOH, n_A=3, pharm_A_idx=0 (OHの酸素)
core     = "CCCC"          # C4リンカー, n_core=4
AnchorB  = "C1CCCCC1"      # シクロヘキシル, n_B=6, pharm_B_offset=0

full_SMILES    = "OC(=O)CCCCC1CCCCC1"
pharm_A_idx    = 0              (OH酸素)
pharm_B_idx    = 3 + 4 + 0 = 7 (シクロヘキシルC1)
```

---

#### 7-5. RDKit距離幾何学 (Distance Geometry) による拘束付き3D埋め込み

距離幾何法（DG法）は分子の内部座標をメトリック行列に変換し、固有値分解により3D座標を求める。RDKitの `rdDistGeom` モジュールが使用する Bounds Matrix（境界行列）を用いて距離制約を課す。

**アルゴリズム**:

```
① GetMoleculeBoundsMatrix(mol_h)
   → 共有結合・角度・vdW 半径から計算したデフォルトバウンズ行列 bm を取得
   　 bm[i][j] = 原子i–j間の上限距離  (i < j の上三角)
   　 bm[j][i] = 原子i–j間の下限距離  (i > j の下三角)

② ファーマコフォア間距離を修正:
   i, j = sorted([pharm_A_idx, pharm_B_idx])
   bm[i][j] = min(d_target + 2.5,  bm[i][j])   # 上限を締める
   bm[j][i] = max(d_target - 2.5,  bm[j][i])   # 下限を広げる

③ EmbedParameters.SetBoundsMat(bm) で制約を組み込み、
   EmbedMultipleConfs(mol_h, numConfs=50) で50配座を生成
   (内部でTriangle Inequality Smoothingを自動適用)

④ MMFFOptimizeMoleculeConfs() で全配座を力場最適化
```

**最良配座の選択**:

制約充足度とMMFFエネルギーの加重スコアで最良配座を選ぶ：

```
score = 0.7 × (|実測距離 - d_target|² / range_dist)
      + 0.3 × (MMFF_energy / range_energy)
```

---

#### 7-6. ドラッグライクネス評価・SDF保存・ドッキング

生成した最良配座について:

1. Lipinski Ro5 (MW≤500, LogP≤5, HBD≤5, HBA≤10) + Veber則 (PSA≤140 Å², RotBonds≤10) でDrugLike判定
2. SDF保存 (`results/bridge/input_bridge_*.sdf`)
3. smina で `--autobox_ligand peptide_ref.pdb` ドッキングを実行

---

#### 7-7. 本PDBでの結果

全ペアで DrugLike 分子が多数生成され、ポケット形状に即した構造を示した。

| ペア | 生成数 | DrugLike数 | 最良ドッキングスコア | 最良SMILES |
|------|-------|-----------|-------------------|------------|
| LYS48–LEU50 | 36 | 36 | **-6.3 kcal/mol** | `O=C(O)CCCCC1CCCCC1` |
| LYS46–TYR49 | 42 | 42 | **-7.2 kcal/mol** | `O=C(O)CC=CCc1ccc2ccccc2c1` |
| ILE21–TYR25 | 48 | 24 | **-7.4 kcal/mol** | `c1ccc2cc(COCCC3CCCCC3)ccc2c1` |

**考察**: POS–ARO型 (LYS46–TYR49) ではナフタレンを用いたアンカーBがTYR49との π–πスタッキングを形成し高スコアを示した。HYD–ARO型 (ILE21–TYR25) でも疎水性コア + 芳香族スタッキングの組み合わせが最良スコアを記録した。

---

### STEP 8: 環状ペプチド vs 低分子 比較

**実装**: `compare_cyclic_peptide.py`

設計した低分子候補が元のペプチドと比べてどの程度の結合力を持つかを定量的に評価するモジュール。

#### 8-1. 環状ペプチド SMILES 自動生成

1文字コードのアミノ酸配列から **head-to-tail 環状ペプチド SMILES** を自動構築する。

```
SMILES 構築ルール (10残基の場合):
  ① 最初の残基 (Gly-1) のバックボーン N に ring closure %11 を付加
  ② 各残基のバックボーン: N-Cα(側鎖)-C(=O) を順に連結
     Pro は特殊: N{r}CCCC{r}C(=O) でピロリジン環を形成
  ③ 最後の残基 (Asp-10) の CO を C(=O)%11 で ring closure

生成例 cyclo-GEVDGWATPD:
  N%11CC(=O)NC(CCC(=O)O)C(=O)NC(C(C)C)C(=O)...N3CCCC3C(=O)NC(CC(=O)O)C(=O)%11
  MW = 1028 Da  |  環サイズ = 30員環マクロラクタム  |  Rings = 4
```

環状ペプチドの特性:
| 指標 | 値 | 注 |
|------|---|----|
| MW | 1028 Da | 30員環マクロラクタム |
| LogP | -4.79 | 親水性 (荷電側鎖多数) |
| HBD / HBA | 14 / 14 | 高 (ペプチド結合由来) |
| RotBonds | 11 | 側鎖のみカウント |

#### 8-2. 3D 配座生成

RDKit の **ETKDG v3** (大環状分子対応) で複数配座を生成し、MMFF94s 最適化後に最低エネルギー配座を選択する。

```python
params = AllChem.ETKDGv3()
params.useMacrocycleTorsions = True   # 大環状回転障壁を考慮
confs = AllChem.EmbedMultipleConfs(mol_h, numConfs=300)
AllChem.MMFFOptimizeMoleculeConfs(mol_h, mmffVariant="MMFF94s")
# → 最低 MMFF エネルギーの配座を SDF 保存
```

#### 8-3. ネイティブポーズ scoring

Chain B の PDB 座標 (`results/docking/peptide_ref.pdb`) を **そのまま** smina に渡し、`--score_only` でスコアを取得する。

```bash
smina --receptor receptor.pdb --ligand peptide_ref.pdb --score_only
→ ネイティブ結合ポーズでの smina Vina スコア (kcal/mol)
```

これは「結晶構造/モデル構造での結合エネルギー」に相当し、低分子のドッキングスコアと直接比較できる。

#### 8-4. 環状ペプチド de novo ドッキング

生成した 3D SDF を smina で通常ドッキングする (低分子と同条件)。

```
exhaustiveness = 16  (大環状は自由度が高いため標準より高めに設定)
autobox_ligand = peptide_ref.pdb  (同じ結合部位)
```

#### 8-5. 比較グラフ

全結果を横棒グラフで可視化 (スコア昇順):

```
色分け:
  ■ 緑 (濃)  : 線形ペプチド ネイティブポーズ (score_only)
  ■ 緑 (淡)  : 環状ペプチド de novo docking
  ■ 青       : ペプチドミメティクス (STEP 3)
  ■ 橙       : ファーマコフォアブリッジ (STEP 7)
  --- 赤破線  : PRODIGY ΔG = -14.7 kcal/mol (参考; 異なるスケール)
```

> **注意**: PRODIGY は protein-peptide 界面の統計モデルに基づく ΔG 予測であり、smina Vina スコアとは定義・スケールが異なります。参考値として表示。

---

### STEP 9: AutoDock CrankPep 環状ペプチドドッキング

**実装**: `dock_cyclic_adcp.py`

smina の代わりに **AutoDock CrankPep (ADCP)** — 環状ペプチド専用ドッキングツール — を用いて
cyclo-GEVDGWATPD をドッキングし、smina 低分子結果と並列比較する。

#### スコアスケールの違い (重要)

| ツール | 手法 | 典型スコア範囲 | 用途 |
|--------|------|--------------|------|
| **ADCP** | Monte Carlo + OpenMM implicit 溶媒 | **-10〜-35 kcal/mol** | 環状ペプチド専用 |
| **smina** | AutoDock Vina スコア関数 | **-5〜-10 kcal/mol** | 低分子 |

> **スコアの直接比較は不可。** 異なるスコアリング関数・溶媒モデルのため数値スケールが異なる。
> LE (Ligand Efficiency = |Affinity| / HAC) も参考値として扱うこと。

#### 9-1. 受容体ターゲットファイル作成 (`--setup-receptor`)

```
receptor.pdb
    │
    ├─ agfr --toPdbqt ──→ receptor_rec.pdbqt
    │
    ├─ obabel -opdbqt ──→ peptide_ref.pdbqt   (結合部位特定用)
    │
    └─ agfr -asv 1.1 ──→ receptor.trg         (AutoSite v1.1 binding pocket)
                          ≈28 MB、初回のみ必要
```

#### 9-2. ADCP 環状ペプチドドッキング

```bash
# 内部実行コマンド
adcp -O -T receptor.trg -s GEVDGWATPD -cyc \
     -N 5 -n 100000 -nmin 5 -nitr 500 -env implicit -dr -reint \
     -o result_gevdgwatpd -w results/adcp_docking/gevdgwatpd/
```

| フラグ | 説明 |
|--------|------|
| `-cyc` | 環状ペプチドモード (N末端−C末端 共有結合) |
| `-N 5` | MC 探索を 5 回実行 |
| `-n 100000` | 各探索で 100,000 ステップ評価 |
| `-nmin 5` | OpenMM エネルギー最小化 5 回 |
| `-nitr 500` | ADCP イテレーション数 |
| `-env implicit` | implicit (暗黙的) 溶媒モデル |
| `-dr -reint` | dry-run 最適化 + re-integration |

出力 (`result_gevdgwatpd_summary.dlg`) から **mode 1 の affinity** を抽出する。

#### 9-3. smina 低分子結果との比較

`Result_Best/summary.csv` (または `results/bridge/bridge_*_results.csv`) から
上位 15 件を読み込み、2 パネルグラフで比較:

- 左パネル: ADCP Affinity + PRODIGY ΔG 参照線
- 右パネル: smina Top15 スコア (LE 色分け)
- LE 比較グラフ (ADCP LE / smina LE 並列)

---

### STEP 10: AutoDockFR 低分子ドッキング

**実装**: `dock_smol_adfr.py`

smina (Vina) の代わりに **AutoDockFR (ADFR)** — ADCP と同じ AutoDock 力場ベースの低分子ドッキングツール — を用いて
Result_Best の候補分子をドッキングし、ADCP 環状ペプチド結果と力場レベルで比較する。

#### スコアスケールの違い (3 メソッド)

| ツール | 対象 | 力場 | 溶媒モデル | 典型スコア範囲 |
|--------|------|------|-----------|--------------|
| **ADCP** | 環状ペプチド | AutoDock 4 | GB/SA implicit | **-10〜-35 kcal/mol** |
| **ADFR** | 低分子 | AutoDock 4 | なし (真空中) | **-5〜-15 kcal/mol** |
| **smina** | 低分子 | Vina 経験的 | なし | **-5〜-10 kcal/mol** |

> **ADFR と ADCP は同じ力場**を使用 (静電相互作用 + 脱溶媒和項 + vdW + HBond)。
> ADCP は GB/SA implicit 溶媒を加えるためスコアが大きくなる傾向がある。
> いずれのメソッドも**スコアの直接比較は不可**。LE も同一メソッド内での比較を推奨。

#### 10-1. 前提条件

```
① adcpsuite 環境のセットアップ (初回のみ)
     bash adcpsuite_micromamba.sh

② receptor.trg の作成 (ADCP と共用 / 初回のみ)
     python dock_cyclic_adcp.py --setup-receptor

③ 低分子候補の選抜 (推奨)
     python collect_best.py
     → Result_Best/summary.csv  ←  ADFR がここから分子を読み込む
```

#### 10-2. SDF → PDBQT 変換

ADFR は PDBQT 形式のリガンドを必要とする。`obabel` (adcpsuite 環境内) で変換する。

```bash
# 内部実行コマンド (フォールバックあり)
obabel -isdf mol.sdf -opdbqt -O mol.pdbqt --partialcharge gasteiger -h
# -h: 極性水素を追加  --partialcharge gasteiger: Gasteiger 部分電荷を計算
# フォールバック 1: -h なしで再試行
# フォールバック 2: prepare_ligand4.py (adcpsuite 付属)
```

#### 10-3. ADFR 低分子ドッキング

```bash
# 内部実行コマンド
adfr -t receptor.trg \
     -l ligand.pdbqt \
     -o result_{mol_name} \
     -N 20 -n 2500000 \
     -w results/adfr_docking/{mol_name}/
```

| フラグ | 説明 |
|--------|------|
| `-t receptor.trg` | ADCP と共用の AutoSite v1.1 ターゲットファイル |
| `-l ligand.pdbqt` | PDBQT 形式リガンド |
| `-N 20` | GA (遺伝的アルゴリズム) 探索回数 |
| `-n 2500000` | 各 GA 実行での評価ステップ数 |
| `--quick` | `-N 5, -n 500000` (高速テスト用) |

出力 (`.dlg` ファイル) から **mode 1 の affinity** を抽出する。

#### 10-4. 生成される比較グラフ (4 種)

| ファイル | 内容 |
|----------|------|
| `adfr_docking_scores.png` | ADFR スコア棒グラフ (LE 色分け: ◎/○/△/✕) |
| `adfr_vs_smina.png` | 同一分子: ADFR (左) vs smina (右) 並列比較 |
| `adfr_vs_adcp_comparison.png` | ADFR 低分子 (右) vs ADCP 環状ペプチド (左) |
| `le_comparison_3methods.png` | LE: ADCP / ADFR / smina 3 メソッド比較バー |

---

## 手法の比較

### 相互作用解析手法の比較

| 手法 | 計算コスト | 動的効果 | 精度 | 実装ファイル |
|------|-----------|----------|------|-------------|
| 距離ベース接触 | 低 | × | 低〜中 | `analyze_interactions.py` |
| ΔSASA | 低 | × | 中 | `analyze_sasa.py` |
| PRODIGY | 低 | × | 中〜高 | `analyze_prodigy.py` |
| smina Vina | 中 | × | 中 | `dock_with_smina.py` |
| smina score_only | 低 | × | 中 | `dock_with_smina.py` |

3手法のスコアが一致する残基ほど、実験的ホットスポットである可能性が高い。

### 低分子生成手法の比較

| 手法 | 設計起点 | 3D形状考慮 | DrugLike担保 | 実装ファイル |
|------|---------|-----------|-------------|-------------|
| ペプチドミメティクス (STEP 3) | ペプチド側鎖フラグメント | × | Ro5チェックのみ | `design_small_molecule.py` |
| ポケット相補設計 + BRICS (STEP 3b) | タンパク質ポケット残基タイプ | △ (テンプレート) | Ro5 + BRICSフィルタ | `generate_pocket_molecules.py` |
| **ファーマコフォアブリッジ (STEP 7)** | **ポケット残基間Cβ距離** | **◎ (DG法拘束)** | Ro5 + Veber | `pharmacophore_bridge.py` |

STEP 7 の距離幾何学アプローチは、ポケット形状を定量的に反映した最初の手法であり、生成された3D構造がより実際の結合部位に適合する。

---

## 使い方

### 必要なもの

| ソフトウェア | バージョン | 備考 |
|---|---|---|
| macOS | 12 Monterey 以降 | `smina.osx.12` の動作要件 |
| Python | 3.9 以上 | 3.14 で動作確認済み |
| micromamba | 最新 | ADCP/ADFR 専用環境用（STEP 9/10 のみ） |

---

### クイックスタート（全体の流れ）

```
① smina 権限設定
        ↓
② Python 仮想環境構築 (venv)
        ↓
③ メインパイプライン実行 (STEP 1〜6)  →  results/ に出力
        ↓
④ STEP 7: ファーマコフォアブリッジ
        ↓
⑤ STEP 8: 環状ペプチド vs 低分子 比較 (smina)
        ↓
⑥ 最終候補の選抜 (collect_best.py)  →  Result_Best/ に出力
        ↓
⑦ micromamba + adcpsuite 環境構築 ← STEP 9/10 に必要（初回のみ）
        ↓
⑧ STEP 9: AutoDock CrankPep 環状ペプチドドッキング (ADCP)
        ↓
⑨ STEP 10: AutoDockFR 低分子ドッキング (ADFR)
```

---

### ① smina バイナリの権限設定（macOS / 初回のみ）

macOS はダウンロードしたバイナリに実行権限がなく、Gatekeeper にもブロックされます。

```bash
# 実行権限を付与
chmod +x smina.osx.12

# Gatekeeper の隔離フラグを解除
xattr -cr smina.osx.12
```

> **確認方法**: `./smina.osx.12 --version` でバージョン文字列が表示されれば OK。
> それでもブロックされる場合は [システム設定 → プライバシーとセキュリティ] で許可してください。

---

### ② Python 仮想環境の構築（初回のみ）

```bash
# 仮想環境を作成（プロジェクト直下に venv/ が生成される）
python3 -m venv venv

# パッケージをインストール
venv/bin/pip install --upgrade pip
venv/bin/pip install biopython numpy scipy matplotlib rdkit
```

> **rdkit が pip でインストールできない場合**（Python バージョン非対応など）:
> ```bash
> # conda または micromamba を使う
> conda create -n peptide_pipeline python=3.11 rdkit -c conda-forge -y
> conda activate peptide_pipeline
> pip install biopython numpy scipy matplotlib
> # 以降のコマンドの venv/bin/python を python に読み替えてください
> ```

---

### ③ メインパイプラインの実行（STEP 1〜6）

```bash
# フル実行（相互作用解析 → 低分子設計 → smina ドッキング）
venv/bin/python pipeline.py Protein_Peptide.pdb
```

**所要時間**: 15〜30 分（デフォルト設定）
**出力先**: `results/`

```bash
# 設計に使う残基を増やす（候補数が増える）
venv/bin/python pipeline.py --top-residues 4

# ドッキング精度を上げる（時間が増える）
venv/bin/python pipeline.py --exhaustiveness 16 --num-modes 9

# 設計だけ行い、ドッキングをスキップ
venv/bin/python pipeline.py --skip-docking

# ドッキングのみ単体実行（既存 SDF を再利用する場合）
venv/bin/python dock_with_smina.py Protein_Peptide.pdb \
    --ligands results/candidate_ligands.sdf \
    --output-dir results/docking
```

---

### ④ 最終候補の選抜

STEP 6 が終わったら上位候補を `Result_Best/` にまとめます。

```bash
venv/bin/python collect_best.py
# → Result_Best/summary.csv  (STEP 10 の ADFR ドッキングで読み込む)
# → Result_Best/01_score*.sdf 〜 15_score*.sdf
```

---

### ⑤ ADCP/ADFR 環境構築（初回のみ・STEP 9/10 を使う場合）

STEP 9 (ADCP) と STEP 10 (ADFR) には Python 3.7 の専用環境 **adcpsuite** が必要です。
メインパイプラインの `venv` とは**別の環境**で、`micromamba` で管理します。

#### 5-1. micromamba のインストール

```bash
# Homebrew を使う場合（推奨）
brew install micromamba

# または公式インストーラ
curl -Ls https://micro.mamba.pm/install.sh | bash
# → ~/.zshrc または ~/.bashrc に初期化コードが追記される
# → ターミナルを再起動するか source ~/.zshrc を実行
```

> **確認方法**: `micromamba --version` でバージョンが表示されれば OK。

#### 5-2. adcpsuite 環境のセットアップ

```bash
bash adcpsuite_micromamba.sh
# → micromamba 環境 adcpsuite (Python 3.7) を作成
# → ADCP v0.0.25 / ADFR v0.0.22 / AGFR / OpenMM / obabel をインストール
# 所要時間: 5〜20 分（ダウンロード込み）
```

#### 5-3. 受容体ターゲットファイルの作成（初回のみ）

ADCP と ADFR は共通の `.trg` ファイルを使います。一度だけ作成すれば以降は不要です。

```bash
python dock_cyclic_adcp.py --setup-receptor
# → results/adcp_docking/receptor.trg を生成 (≈28 MB)
# 所要時間: 2〜5 分
```

---

### STEP 7: ファーマコフォアブリッジ分子生成

```bash
# 既知ポケットペアを一覧表示
venv/bin/python pharmacophore_bridge.py --list-pairs

# LYS48–LEU50 間を橋渡しする分子を生成 + ドッキング
venv/bin/python pharmacophore_bridge.py --point1 LYS48 --point2 LEU50

# LYS46–TYR49 をターゲット（ドッキングなし、高速確認）
venv/bin/python pharmacophore_bridge.py --point1 LYS46 --point2 TYR49 --no-dock

# 任意残基ペアを指定（ARG28–TYR49, 9.1 Å）
venv/bin/python pharmacophore_bridge.py --point1 ARG28 --point2 TYR49

# コンフォマー数を増やして探索を強化
venv/bin/python pharmacophore_bridge.py \
    --point1 ILE21 --point2 TYR25 \
    --n-confs 100 --exhaustiveness 16
```

#### STEP 7 オプション一覧

| オプション | デフォルト | 説明 |
|---|---|---|
| `--point1` | — | 残基1 (例: `LYS48`) |
| `--point2` | — | 残基2 (例: `LEU50`) |
| `--chain` | `A` | タンパク質チェーンID |
| `--no-dock` | — | ドッキングをスキップ |
| `--n-confs` | `50` | 生成コンフォマー数 |
| `--tolerance` | `2.5` | 距離制約の許容範囲 (Å) |
| `--exhaustiveness` | `8` | smina 探索徹底度 |
| `--list-pairs` | — | 既知ポケットペアを表示して終了 |
| `--all-pairs` | — | 全既知ペアを一括処理 |

### STEP 8: 環状ペプチド vs 低分子 比較

```bash
# 環状ペプチドをドッキングして既存低分子と比較 (推奨)
venv/bin/python compare_cyclic_peptide.py Protein_Peptide.pdb

# 以前のドッキング結果を使って比較グラフのみ再生成 (高速)
venv/bin/python compare_cyclic_peptide.py --no-dock

# 配座生成数を増やして精度向上 (時間増)
venv/bin/python compare_cyclic_peptide.py --n-confs 500 --exhaustiveness 24

# 異なるペプチド配列で比較
venv/bin/python compare_cyclic_peptide.py --sequence GEVDGWATPD
```

#### STEP 8 オプション一覧

| オプション | デフォルト | 説明 |
|---|---|---|
| `pdb` | `Protein_Peptide.pdb` | 入力 PDB |
| `--sequence` | `GEVDGWATPD` | ペプチド配列 (1文字コード) |
| `--no-dock` | — | 既存結果からグラフのみ生成 |
| `--n-confs` | `300` | 3D 配座生成数 |
| `--exhaustiveness` | `16` | smina 探索の徹底度 |
| `--results-dir` | `results` | 既存結果ディレクトリ |
| `--output-dir` | `results/comparison` | 出力先ディレクトリ |

---

### STEP 9: AutoDock CrankPep 環状ペプチドドッキング

> **前提**: 上記 ⑤ の adcpsuite 環境構築と receptor.trg の作成が完了していること。

```bash
# 通常実行 (N=5, n=100000 ≈5〜10分)
python dock_cyclic_adcp.py

# クイックテスト (N=1, n=5000 ≈30秒〜1分 / 精度低)
python dock_cyclic_adcp.py --quick

# 既存結果を強制上書きして再実行
python dock_cyclic_adcp.py --force

# 別の配列で実行
python dock_cyclic_adcp.py --sequence ACDEFGHIKL
```

#### STEP 9 オプション一覧

| オプション | デフォルト | 説明 |
|---|---|---|
| `--sequence` | `GEVDGWATPD` | ペプチド配列 (1文字コード) |
| `--setup-receptor` | — | `.trg` ファイルを再作成する |
| `--quick` | — | N=1, n=5000 クイックモード |
| `--n-runs` | `5` | MC 探索回数 |
| `--n-evals` | `100000` | 評価ステップ数 |
| `--timeout` | `3600` | タイムアウト秒数 |
| `--force` | — | 既存結果を上書きして再実行 |

---

### STEP 10: AutoDockFR 低分子ドッキング

```bash
# 通常実行 (Result_Best 上位 15 件を ADFR でドッキング)
python dock_smol_adfr.py

# クイックテスト (N=5, n=500000, ≈1〜2 分/分子 / 精度低)
python dock_smol_adfr.py --quick

# 処理件数を絞る (上位 5 件のみ)
python dock_smol_adfr.py --top-n 5

# 既存結果を上書きして再実行
python dock_smol_adfr.py --force

# GA 探索回数・ステップ数をカスタム指定
python dock_smol_adfr.py --n-runs 30 --n-evals 5000000
```

#### STEP 10 オプション一覧

| オプション | デフォルト | 説明 |
|---|---|---|
| `--top-n` | `15` | ドッキングする分子数 |
| `--quick` | — | N=5, n=500000 クイックモード |
| `--n-runs` | `20` | GA 探索回数 |
| `--n-evals` | `2500000` | 評価ステップ数 |
| `--timeout` | `600` | 1 分子あたりのタイムアウト秒数 |
| `--force` | — | 既存結果を上書きして再実行 |

> **前提**: `bash adcpsuite_micromamba.sh` で adcpsuite 環境を構築し、
> `python dock_cyclic_adcp.py --setup-receptor` で `receptor.trg` を作成済みであること。

---

### オプション一覧

| オプション | デフォルト | 説明 |
|---|---|---|
| `pdb` | `Protein_Peptide.pdb` | 入力PDBファイル |
| `--protein-chain` | `A` | タンパク質チェーンID |
| `--peptide-chain` | `B` | ペプチドチェーンID |
| `--top-residues` | `3` | 低分子設計に使う上位残基数 |
| `--cutoff` | `4.5` | 接触検出の距離カットオフ (Å) |
| `--exhaustiveness` | `8` | sminaの探索の徹底度 |
| `--num-modes` | `9` | sminaが出力するポーズ数 |
| `--skip-docking` | — | ドッキングをスキップ |
| `--skip-rescore` | — | score_only rescoring をスキップ |
| `--skip-sasa` | — | ΔSASA解析をスキップ |
| `--skip-prodigy` | — | PRODIGY予測をスキップ |
| `--output-dir` | `results` | 出力ディレクトリ |

---

## 出力ファイル

```
results/
├── contacts.json              # 全接触データ（距離・座標・タイプ）
├── sasa_analysis.csv          # 残基別ΔSASA（単体/複合体）
├── sasa_comparison.png        # 埋没面積グラフ
├── prodigy_result.json        # PRODIGY ΔG・Kd予測結果
├── prodigy_contacts.png       # 界面接触タイプ円グラフ
├── pharmacophore.csv          # ファーマコフォア特徴点（x,y,z座標）
├── pharmacophore.pml          # PyMOL用ファーマコフォア描画スクリプト
├── candidate_ligands.sdf      # 設計した低分子候補（3D配座つき）
├── report.json                # パイプライン全体のサマリー（PRODIGY/SASA含む）
├── residue_scores.png         # 残基重要度グラフ
├── interaction_map.png        # 相互作用ヒートマップ
├── pharmacophore_3d.png       # ファーマコフォア3D図
└── docking/
    ├── receptor.pdb           # 受容体（Chain A）
    ├── peptide_ref.pdb        # 参照ペプチド（Chain B, autobox用）
    ├── docked_*.sdf           # ドッキングポーズ（PyMOLで可視化可）
    ├── log_*.txt              # sminaログ（全ポーズのスコア）
    ├── docking_results.csv    # ドッキングスコアランキング
    ├── docking_scores.png     # スコア比較グラフ
    ├── rescore_terms.csv      # score_only スコア項分解（全ポーズ）
    ├── score_decomposition.png# スコア項分解グラフ（gauss/HYD/HB/REP）
    └── rescoring/
        ├── pose_*.sdf         # ポーズ別個別SDF
        └── rescore_*.txt      # score_only ログ
```

```
results/bridge/                           # STEP 7 (ファーマコフォアブリッジ) 出力
├── bridge_LYS48_LEU50_results.csv        # 生成分子のスコア・性質一覧
├── bridge_LYS46_TYR49_results.csv
├── bridge_ILE21_TYR25_results.csv
├── bridge_*_analysis.png                 # ドッキングスコア + 距離充足グラフ
├── geometry_*.png                        # Cβ–Cβ距離の3D俯瞰図
├── input_bridge_*.sdf                    # ドッキング入力SDF（3D配座つき）
└── docked_bridge_*.sdf                   # ドッキングポーズ
```

CSV の主要カラム:

| カラム | 内容 |
|--------|------|
| `smiles` | 正規化SMILES |
| `res1`, `res2` | 対象残基ペア |
| `anchor_a`, `anchor_b` | 使用アンカー基名 |
| `linker` | 使用リンカー名 |
| `target_dist` | 目標Cβ–Cβ距離 (Å) |
| `actual_dist` | 最良配座での実測ファーマコフォア間距離 (Å) |
| `dock_score` | sminaドッキングスコア (kcal/mol) |
| `MW`, `LogP`, `HBD`, `HBA`, `PSA`, `RotBonds` | 物理化学的性質 |
| `DrugLike` | Ro5 + Veber 充足フラグ |

```
results/comparison/                           # STEP 8 (smina 環状ペプチド比較) 出力
├── cyclic_peptide_GEVDGWATPD.sdf             # 環状ペプチド 3D 配座 (ドッキング入力)
├── docked_cyclo_GEVDGWATPD.sdf              # 環状ペプチド ドッキングポーズ
├── docked_cyclo_GEVDGWATPD.log              # smina ドッキングログ
├── native_peptide_score.log                 # 線形ペプチド score_only ログ
├── binding_comparison.png                   # 結合親和性比較グラフ
├── le_comparison.png                        # LE 比較グラフ
└── comparison_summary.json                  # 比較結果サマリー (JSON)
```

```
results/adcp_docking/                         # STEP 9 (AutoDock CrankPep) 出力
├── receptor.pdb                              # 受容体 PDB (コピー)
├── receptor_rec.pdbqt                        # 受容体 PDBQT (agfr 変換)
├── peptide_ref.pdbqt                         # ペプチド PDBQT (obabel 変換)
├── receptor.trg                              # ADCP ターゲットファイル (≈28 MB, agfr 作成)
├── gevdgwatpd/                               # ADCP ドッキング結果
│   ├── result_gevdgwatpd_summary.dlg         # ドッキングサマリー (affinity ここから)
│   └── result_gevdgwatpd_out.pdb             # 最良ポーズ 3D 構造 (PyMOL 可視化可)
├── adcp_vs_smina_comparison.png              # 2 パネル比較グラフ (ADCP vs smina)
├── adcp_vs_smina_le_comparison.png           # LE 比較グラフ
└── adcp_comparison_report.json               # 結果サマリー (JSON)
```

```
# adcp_comparison_report.json の主要フィールド
{
  "adcp_cyclic_peptide": {
    "sequence": "GEVDGWATPD",
    "affinity_kcal_mol": -XX.X,   # ADCP スコア (典型: -10〜-35)
    "hac": 73,
    "le_adcp": 0.XXXX             # LE = |affinity| / HAC (参考値)
  },
  "smina_top_small_molecules": [...],
  "note": "ADCP と smina は異なる scoring function のため直接比較不可"
}
```

```
results/adfr_docking/                         # STEP 10 (AutoDockFR) 出力
├── ligands_pdbqt/                            # SDF → PDBQT 変換済みファイル
│   ├── bridge_LYS46_TYR49_*.pdbqt
│   └── ...
├── {mol_name}/                               # 各分子の ADFR 出力ディレクトリ
│   ├── result_{name}_summary.dlg             # ドッキングサマリー (affinity ここから)
│   └── result_{name}_out.pdbqt              # 最良ポーズ 3D 構造 (PDBQT)
├── adfr_docking_scores.png                   # ADFR スコア棒グラフ (LE 色分け)
├── adfr_vs_smina.png                         # 同一分子: ADFR vs smina 並列比較
├── adfr_vs_adcp_comparison.png               # ADFR 低分子 vs ADCP 環状ペプチド比較
├── le_comparison_3methods.png               # LE: ADCP / ADFR / smina 3 メソッド比較
├── adfr_docking_results.csv                  # 全結果 CSV (ADFR+smina スコア/LE)
└── adfr_comparison_report.json              # 結果サマリー JSON
```

```
# adfr_comparison_report.json の主要フィールド
{
  "adfr_small_molecules": [
    {
      "rank": 1,
      "name": "bridge_LYS46_TYR49_...",
      "adfr_score": -X.XXX,       # ADFR スコア (典型: -5〜-15 kcal/mol)
      "smina_score": -X.XXX,      # smina スコア (参照用)
      "hac": XX,
      "le_adfr": 0.XXXX,          # LE = |ADFR score| / HAC
      "smina_le": 0.XXXX,
      "pair": "LYS46-TYR49"
    }, ...
  ],
  "adcp_cyclic_peptide": { ... }, # STEP 9 と共用 (存在する場合)
  "scoring_info": {
    "ADFR":  { "force_field": "AutoDock 4", "solvent": "なし", "typical_range": "-5〜-15" },
    "ADCP":  { "force_field": "AutoDock 4 + GB/SA", "typical_range": "-10〜-35" },
    "smina": { "force_field": "Vina 経験的", "typical_range": "-5〜-10" },
    "comparison_note": "異なるメソッド間の直接スコア比較は不可"
  }
}
```

### PyMOLでの可視化

```python
# PyMOL上で実行: 基本ドッキング結果
load results/docking/receptor.pdb
load results/docking/peptide_ref.pdb
load results/docking/docked_2_reduced_peptide_TRP_.sdf
run  results/pharmacophore.pml

# STEP 7 ブリッジ分子の可視化
load results/docking/receptor.pdb
load results/bridge/docked_bridge_LYS46_TYR49_carboxylic_acid_diene_naphthalene.sdf

# STEP 9 ADCP 環状ペプチド最良ポーズの可視化
load results/docking/receptor.pdb
load results/adcp_docking/gevdgwatpd/result_gevdgwatpd_out.pdb
show cartoon, result_gevdgwatpd_out
color cyan, result_gevdgwatpd_out

# STEP 10 ADFR 低分子ドッキング最良ポーズの可視化
load results/docking/receptor.pdb
load results/adfr_docking/bridge_LYS46_TYR49_.../result_bridge_LYS46_TYR49_..._out.pdbqt
show sticks, result_bridge_LYS46_TYR49_..._out
color magenta, result_bridge_LYS46_TYR49_..._out
```

---

## 最終選抜結果 (Result_Best)

全ドッキング結果から **スコア上位 + Ligand Efficiency 上位** を選抜した最終候補。
`collect_best.py` を実行すると `Result_Best/` に自動生成される。

```bash
python collect_best.py
# → Result_Best/ に上位 SDF + summary.csv + README.md を生成
```

### 選抜基準

| 基準 | 内容 |
|------|------|
| DrugLike | Lipinski Ro5 + Veber 充足 |
| スコアカットオフ | ≤ -5.0 kcal/mol |
| スコア上位 | Top 10 (より負 = 強い結合) |
| LE 上位 | Top 5 (重複除去後に追加) |

### LE (Ligand Efficiency) の見方

```
LE = |smina Affinity| / 重原子数 (HAC)

◎ LE ≥ 0.4  : 優秀 (フラグメント〜低分子の理想域)
○ LE ≥ 0.3  : 良好 (低分子医薬品の目安)
△ LE ≥ 0.2  : 許容範囲
✕ LE < 0.2  : 非効率 (大きすぎ or 弱すぎ)
```

### 選抜分子一覧 (本PDB 結果)

| 順位 | スコア (kcal/mol) | HAC | LE | LE評価 | 対象残基ペア | 分子名 |
|------|-----------------|-----|-----|--------|------------|--------|
| 1 | -7.200 | 26 | 0.2769 | △ 許容 | TRP6-ALA7-ASP4 | reduced_peptide_TRP_ALA_PRO |
| 2 | -7.200 | 17 | 0.4235 | ◎ 優秀 | LYS46-TYR49 | bridge_LYS46_TYR49_carboxylic_acid_diene_naphthalene |
| 3 | -7.200 | 18 | 0.4000 | ◎ 優秀 | LYS46-TYR49 | bridge_LYS46_TYR49_sulfonamide_diene_naphthalene |
| 4 | -7.100 | 20 | 0.3550 | ○ 良好 | LYS46-TYR49 | bridge_LYS46_TYR49_sulfonamide_amide_lnk_naphthalene |
| 5 | -7.000 | 17 | 0.4118 | ◎ 優秀 | LYS46-TYR49 | bridge_LYS46_TYR49_carboxylic_acid_ether3_naphthalene |
| 6 | -7.000 | 19 | 0.3684 | ○ 良好 | LYS46-TYR49 | bridge_LYS46_TYR49_carboxylic_acid_amide_lnk_naphthalene |
| 7 | -7.000 | 18 | 0.3889 | ○ 良好 | LYS46-TYR49 | bridge_LYS46_TYR49_sulfonamide_ether3_naphthalene |
| 8 | -6.900 | 17 | 0.4059 | ◎ 優秀 | LYS46-TYR49 | bridge_LYS46_TYR49_carboxylic_acid_butylene_naphthalene |
| 9 | -6.700 | 18 | 0.3722 | ○ 良好 | LYS46-TYR49 | bridge_LYS46_TYR49_carboxylic_acid_pentylene_naphthalene |
| 10 | -6.700 | 19 | 0.3526 | ○ 良好 | LYS46-TYR49 | bridge_LYS46_TYR49_carboxylic_acid_hexylene_naphthalene |
| 11 | -6.300 | 13 | 0.4846 | ◎ 優秀 | LYS48-LEU50 | bridge_LYS48_LEU50_carboxylic_acid_butylene_cyclohexyl |
| 12 | -6.300 | 13 | 0.4846 | ◎ 優秀 | LYS48-LEU50 | bridge_LYS48_LEU50_carboxylic_acid_butylene_phenyl |
| 13 | -5.600 | 10 | 0.5600 | ◎ 優秀 | LYS48-LEU50 | bridge_LYS48_LEU50_carboxylic_acid_butylene_isopropyl |
| 14 | -5.500 | 10 | 0.5500 | ◎ 優秀 | TRP6-ALA7-ASP4 | fragment_TRP |
| 15 | -5.400 | 11 | 0.4909 | ◎ 優秀 | LYS48-LEU50 | bridge_LYS48_LEU50_carboxylic_acid_pentylene_isopropyl |

> **最高LE**: `bridge_LYS48_LEU50_carboxylic_acid_butylene_isopropyl` (LE=0.560, score=-5.6 kcal/mol)
> **最強結合**: `bridge_LYS46_TYR49_carboxylic_acid_diene_naphthalene` (score=-7.2 kcal/mol, LE=0.4235)

### PyMOL での可視化

```python
load results/docking/receptor.pdb
load Result_Best/01_score*.sdf
# 受容体を表面表示
show surface, receptor
set transparency, 0.5, receptor
```

### 注意事項

- smina (AutoDock Vina) スコアは**同一受容体での比較**に有効
- ペプチド vs 低分子の直接比較は LE を参照すること
- 実験的検証 (IC50, SPR, ITC 等) で確認が必要

---

## 依存関係

### メインパイプライン環境 (venv)

| パッケージ | バージョン | 用途 |
|---|---|---|
| Python | 3.9 以上 (3.14 確認済み) | 実行環境 |
| biopython | 1.86 | PDB解析・NeighborSearch |
| rdkit | 2025.9.5 | 分子操作・3D配座生成 |
| numpy | 2.4.2 | 座標計算 |
| scipy | 1.17.0 | 数値計算 |
| matplotlib | 3.10.8 | グラフ生成 |
| smina | dc3dfab+ | ドッキング（AutoDock Vina 1.1.2ベース、macOS バイナリ同梱） |

```bash
# インストールコマンド
venv/bin/pip install biopython numpy scipy matplotlib rdkit
```

### ADCP/ADFR 環境 (adcpsuite / micromamba)

| パッケージ | バージョン | 用途 |
|---|---|---|
| Python | 3.7 | adcpsuite 専用環境 |
| adcp | 0.0.25 | AutoDock CrankPep（環状ペプチドドッキング） |
| adfr | 0.0.22 | AutoDockFR（低分子ドッキング） |
| agfr | — | 受容体 PDBQT / .trg ファイル作成 |
| openmm | 7.6.0 | GB/SA implicit 溶媒計算（ADCP 内部使用） |
| openbabel | 2.4.1 | SDF → PDBQT 変換 |

```bash
# インストールコマンド（micromamba が必要）
bash adcpsuite_micromamba.sh
```
