# Peptide to Small Molecule Conversion Platform

ペプチド-タンパク質複合体の**結合相互作用を保持**しながら、ペプチドを低分子リガンドへ変換する計算パイプラインです。

入力は PDB ファイル 1 つだけ。相互作用解析 → 低分子設計 → ドッキング検証 → 合成ルート提案までを自動で行います。

---

## 目次

1. [クイックスタート](#クイックスタート)
2. [環境構築](#環境構築)
3. [使い方](#使い方)
4. [パイプラインの全体像](#パイプラインの全体像)
5. [各ステップの解説](#各ステップの解説)
6. [出力ファイル一覧](#出力ファイル一覧)
7. [最終選抜結果 (Result_Best)](#最終選抜結果-result_best)
8. [ファイル構成](#ファイル構成)
9. [トラブルシューティング](#トラブルシューティング)
10. [ライセンス](#ライセンス)

---

## クイックスタート

最短手順で動かしたい方向けのガイドです。

```bash
# 1. リポジトリをクローン
git clone https://github.com/Barashin/Peptide_to_small_molecule.git
cd Peptide_to_small_molecule

# 2. conda 環境を作成 (Python 3.11 + RDKit)
conda create -n peptide_pipeline python=3.11 rdkit -c conda-forge -y
conda activate peptide_pipeline

# 3. Python パッケージをインストール
pip install biopython numpy scipy matplotlib pillow

# 4. smina (ドッキングエンジン) をインストール
conda install -c conda-forge smina -y

# 5. パイプラインを実行 (同梱のサンプル PDB を使用)
python pipeline.py Protein_Peptide.pdb

# 6. 結果の確認
#    results/          → 全解析結果
#    Result_Best/      → 上位候補の SDF + 合成スキーム図
```

> **所要時間の目安**: サンプル PDB で約 10〜30 分 (マシンスペックにより変動)

---

## 環境構築

### 前提条件

- **OS**: Linux (推奨) / macOS / WSL2
- **conda** または **mamba**: [Miniconda](https://docs.conda.io/en/latest/miniconda.html) か [Miniforge](https://github.com/conda-forge/miniforge) をインストール済みであること

### Step 1: conda 環境の作成

```bash
# 新しい環境を作成 (名前は自由に変更可)
conda create -n peptide_pipeline python=3.11 rdkit -c conda-forge -y

# 環境を有効化 (以降のコマンドはすべてこの環境内で実行)
conda activate peptide_pipeline
```

> **なぜ conda ?**: RDKit は pip だけではインストールが難しいため、conda-forge を使います。

### Step 2: Python パッケージのインストール

```bash
pip install biopython numpy scipy matplotlib pillow
```

| パッケージ | バージョン目安 | 用途 |
|-----------|-------------|------|
| biopython | ≥1.80 | PDB ファイル解析・相互作用検出 |
| rdkit | ≥2023 | 分子操作・3D 配座生成・BRICS 分解 |
| numpy | ≥1.24 | 座標計算 |
| scipy | ≥1.10 | 数値計算 |
| matplotlib | ≥3.7 | グラフ描画 |
| pillow | ≥9.0 | 合成スキーム画像生成 |

### Step 3: smina のインストール

smina は分子ドッキングエンジンです。以下のいずれかの方法でインストールしてください。

```bash
# 方法 A: conda-forge からインストール (推奨)
conda install -c conda-forge smina -y

# 方法 B: 静的バイナリをダウンロード (Linux)
wget https://sourceforge.net/projects/smina/files/smina.static/download -O smina.static
chmod +x smina.static
# パスを通すか、実行時に --smina-path で指定
```

動作確認:
```bash
smina --version
# "smina ..." と表示されれば OK
```

### Step 4 (オプション): AiZynthFinder のインストール

AI ベースの逆合成ルート提案機能を使う場合にインストールします。**なくても基本パイプラインは動きます**。

```bash
# AiZynthFinder 本体をインストール
pip install aizynthfinder[all]

# pytables を conda-forge から入れ直す (blosc2 ライブラリの問題を回避)
conda install -c conda-forge pytables -y

# 学習済みモデル (約 750 MB) をダウンロード
download_public_data aizynthfinder_data
```

動作確認:
```bash
python -c "from aizynthfinder.aizynthfinder import AiZynthFinder; print('OK')"
```

> **注意**: `download_public_data` で取得するモデルデータは約 750 MB あります。ダウンロードには時間がかかる場合があります。

### Step 5 (オプション): ADCP/ADFR 環境

環状ペプチドドッキング (AutoDock CrankPep) と AutoDockFR 低分子ドッキングを使う場合のみ必要です。

```bash
# 同梱スクリプトで micromamba 環境をセットアップ
bash adcpsuite_micromamba.sh
```

> これは `peptide_pipeline` とは**別の Python 3.7 環境**を作ります。パイプラインから自動で呼び出されます。

---

## 使い方

### 基本的な実行

```bash
conda activate peptide_pipeline
python pipeline.py <PDBファイル>
```

同梱のサンプルで試す場合:
```bash
python pipeline.py Protein_Peptide.pdb
```

### よく使うオプション

```bash
# チェーンを指定 (デフォルト: A=タンパク質, B=ペプチド)
python pipeline.py input.pdb --protein-chain A --peptide-chain B

# 出力先を変更
python pipeline.py input.pdb --output-dir my_results

# ドッキング精度を上げる (時間がかかる)
python pipeline.py input.pdb --exhaustiveness 16

# AiZynthFinder で AI 逆合成を有効化
python pipeline.py input.pdb --use-aizynthfinder

# ドッキングをスキップ (相互作用解析と低分子設計のみ)
python pipeline.py input.pdb --skip-docking
```

### 全オプション一覧

```
python pipeline.py --help

引数:
  pdb                    入力 PDB ファイル (デフォルト: Protein_Peptide.pdb)

オプション:
  --protein-chain        タンパク質チェーン ID (デフォルト: A)
  --peptide-chain        ペプチドチェーン ID (デフォルト: B)
  --top-residues N       設計に使う上位残基数 (デフォルト: 3)
  --cutoff DIST          接触検出カットオフ距離 [Å] (デフォルト: 4.5)
  --output-dir DIR       出力ディレクトリ (デフォルト: results)
  --exhaustiveness N     smina 探索精度 (デフォルト: 8)
  --num-modes N          ドッキングポーズ数 (デフォルト: 9)
  --sa-threshold SCORE   SA Score 合成容易性閾値 (デフォルト: 6.0)
  --skip-docking         ドッキングをスキップ
  --skip-rescore         リスコアリングをスキップ
  --skip-sasa            ΔSASA 解析をスキップ
  --skip-prodigy         PRODIGY 解析をスキップ
  --skip-synth-filter    合成容易性フィルタをスキップ
  --skip-rigid-scaffold  剛直スキャフォールド設計をスキップ
  --use-aizynthfinder    AiZynthFinder AI 逆合成を有効化
  --aizynthfinder-config AiZynthFinder config.yml のパス
```

### 最終選抜の実行

パイプライン実行後に、全結果から上位候補を選抜して合成スキーム図を生成します:

```bash
python collect_best.py
# → Result_Best/ に Top 5 の SDF + 合成スキーム PNG + HTML レポート を出力
```

### 自分の PDB を使う場合

1. PDB ファイルに**タンパク質とペプチドが別チェーン**で含まれていることを確認
2. チェーン ID を `--protein-chain` と `--peptide-chain` で指定

```bash
# 例: タンパク質が Chain C、ペプチドが Chain D の場合
python pipeline.py my_complex.pdb --protein-chain C --peptide-chain D
```

---

## パイプラインの全体像

```
入力: Protein_Peptide.pdb (タンパク質-ペプチド複合体)
        │
        ▼
┌─────────────────────────────────────────────────────────┐
│  STEP 1   相互作用解析                                   │
│           距離ベース接触検出 → H結合/疎水性/静電 分類      │
│           ΔSASA 解析 → 界面埋没面積の定量                 │
│           PRODIGY → 結合親和性 (ΔG) 予測                  │
│           → 残基重要度ランキング (上位3残基を自動選択)     │
├─────────────────────────────────────────────────────────┤
│  STEP 2   ファーマコフォア抽出                           │
│           重要残基の側鎖 → 薬理活性特徴点 (HBA/HBD/HYD等) │
├─────────────────────────────────────────────────────────┤
│  STEP 3   低分子設計                                     │
│           ペプチドミメティクス (4戦略で候補分子を自動生成)  │
│           合成容易性評価 (SA Score / QED / PAINS / BRENK)│
│           逆合成解析 (BRICS / RECAP / AiZynthFinder)     │
├─────────────────────────────────────────────────────────┤
│  STEP 4   ドッキング検証                                 │
│           smina (AutoDock Vina 系) でタンパク質ポケットに │
│           候補分子をドッキング → 結合スコアでランキング    │
├─────────────────────────────────────────────────────────┤
│  STEP 5   結果集約・選抜                                 │
│           スコア上位 + Ligand Efficiency 上位で Top 5 選抜│
│           合成スキーム図 (AiZynthFinder) + HTML レポート  │
└─────────────────────────────────────────────────────────┘
        │
        ▼
出力: results/          (全解析結果)
      Result_Best/      (上位候補 SDF + 合成スキーム図)
```

---

## 各ステップの解説

### STEP 1: 相互作用解析

ペプチドのどの残基がタンパク質との結合に最も重要かを特定します。

#### 1-1. 距離ベース相互作用解析 (`analyze_interactions.py`)

ペプチドの各原子から **4.5 A** 以内のタンパク質原子を検出し、接触タイプ (H結合 / 疎水性 / 静電 / vdW) に分類します。各ペプチド残基に対して、接触タイプの加重和でスコアを算出します。

```
重み:  H結合 = 3.0 / 静電 = 2.0 / 疎水性 = 1.0 / vdW = 0.5
```

#### 1-2. ΔSASA 解析 (`analyze_sasa.py`)

複合体形成による**溶媒露出面積の変化量** (ΔSASA) を各残基について計算します。
ΔSASA が大きい残基ほど界面に深く埋没しており、結合への寄与が大きいと推定されます。

#### 1-3. PRODIGY 結合親和性予測 (`analyze_prodigy.py`)

Vangone & Bonvin (eLife, 2015) の手法を再実装し、残基間接触タイプのカウントから結合自由エネルギー (ΔG) と解離定数 (Kd) を予測します。

#### スコア統合

距離スコアとΔSASA を正規化して加重平均 (0.6:0.4) で統合し、残基重要度ランキングを生成します。上位残基を STEP 2 以降の低分子設計に使用します。

---

### STEP 2: ファーマコフォア抽出 (`extract_pharmacophore.py`)

重要残基の側鎖から**薬理活性特徴点** (ファーマコフォア) を抽出します。

| 特徴タイプ | 略称 | 検出対象 |
|-----------|------|---------|
| 水素結合受容体 | HBA | O, OD1, OE1 等 |
| 水素結合供与体 | HBD | N, ND1, NZ 等 |
| 疎水性 | HYD | 側鎖 C 原子 |
| 芳香族 | ARO | 芳香環の重心 |
| 正電荷 | POS | LYS の NZ, ARG の NH1/NH2 |
| 負電荷 | NEG | ASP の OD1/OD2, GLU の OE1/OE2 |

---

### STEP 3: 低分子設計 (`design_small_molecule.py`)

重要残基の側鎖構造をもとに、4 つの戦略で候補分子を自動生成します。

| 戦略 | 概要 | 特徴 |
|------|------|------|
| フラグメント単体 | 最重要残基の側鎖をそのまま使用 | 最小サイズ、高 LE |
| 距離適応リンカー連結 | Cβ-Cβ 距離に合わせたリンカーで側鎖を結合 | 柔軟性あり |
| レデュースドペプチド | アミド骨格を残しつつ非ペプチド化 | ペプチドに近い |
| スキャフォールドテンプレート | ベンゼン環等の剛直骨格に側鎖を配置 | Drug-like |

#### 合成容易性評価 (`synthesizability.py`)

| 指標 | 範囲 | 合格基準 |
|------|------|---------|
| SA Score | 1〜10 (低い=容易) | ≤ 6.0 |
| QED | 0〜1 (高い=良好) | ≥ 0.2 |
| PAINS | OK / NG | OK (偽陽性構造なし) |
| BRENK | OK / NG | OK (有害構造なし) |

#### 逆合成解析 (`retrosynthesis.py`)

候補分子の**合成ルートの実現可能性**を推定します。

| 手法 | 説明 |
|------|------|
| BRICS | 合成的に意味のある 16 種の結合で切断 → フラグメントの市販品入手性を評価 |
| RECAP | メディシナルケミストリーの一般的反応で切断 (BRICS と相補的) |
| AiZynthFinder (オプション) | AI ベース逆合成。MCTS + ニューラルネット (USPTO 40K テンプレート) で化学的に正しい合成ルートを提案 |

---

### STEP 3 追加: ファーマコフォアブリッジ (`pharmacophore_bridge.py`)

ポケット内の**タンパク質側 2 残基を同時に橋渡し**できる低分子を設計します。
Cβ-Cβ 距離を拘束条件として、RDKit の距離幾何学 (DG 法) で 3D 配座を生成します。

### STEP 3 追加: 剛直スキャフォールド設計 (`rigid_scaffold_design.py`)

ベンゼン環やナフタレンなどの**剛直な芳香族骨格**を使って、残基間距離にマッチする分子を設計します。conformance rate (距離適合率) で評価し、回転可能結合が少ない Drug-like な候補を優先します。

---

### STEP 4: ドッキング検証 (`dock_with_smina.py`)

smina (AutoDock Vina 系エンジン) で候補分子をタンパク質ポケットにドッキングし、結合親和性を評価します。

- 参照ペプチドの座標からドッキングボックスを自動定義
- 各候補分子を個別にドッキング → アフィニティスコア (kcal/mol) で ランキング
- リスコアリング (`--score_only`) でスコア項分解 (H結合/疎水性/立体反発) も実施

---

### STEP 5: 結果集約・最終選抜 (`collect_best.py`)

全結果を集約し、**ドッキングスコア上位 + Ligand Efficiency 上位**で Top 5 を選抜します。

各候補に対して:
- 逆合成スキーム図 (PNG) を生成
- AiZynthFinder のルートがあれば、化学的に正しい前向き合成シーケンスを描画
- HTML レポートを生成 (構造式 + 逆合成結果 + 合成容易性評価を一覧表示)

---

## 出力ファイル一覧

### メイン出力 (`results/`)

```
results/
├── contacts.json               # 相互作用データ (JSON)
├── residue_scores.png          # 残基重要度スコアの棒グラフ
├── interaction_map.png         # 接触数ヒートマップ
├── pharmacophore.csv           # ファーマコフォア特徴点の 3D 座標
├── pharmacophore.pml           # PyMOL 可視化スクリプト
├── pharmacophore_3d.png        # ファーマコフォア 3D 散布図
├── candidate_ligands.sdf       # 候補分子の SDF
├── docking/                    # ドッキング結果
│   ├── receptor.pdb            #   受容体構造
│   ├── peptide_ref.pdb         #   参照ペプチド
│   ├── docked_*.sdf            #   ドッキングポーズ
│   ├── docking_results.csv     #   スコア一覧
│   └── docking_scores.png      #   スコアグラフ
├── bridge/                     # ファーマコフォアブリッジ結果
├── rigid_scaffold/             # 剛直スキャフォールド結果
├── comparison/                 # 環状ペプチド比較結果
├── adcp_docking/               # AutoDock CrankPep 結果
└── adfr_docking/               # AutoDockFR 結果
```

### 最終選抜 (`Result_Best/`)

```
Result_Best/
├── 01_score*_<name>.sdf        # 選抜分子のドッキングポーズ (PyMOL で可視化可)
├── 02_score*_<name>.sdf
├── ...
├── summary.csv                 # 全選抜分子の詳細データ
├── README.md                   # 選抜基準・分子一覧
└── retrosynthesis/             # 逆合成解析結果
    ├── 01_retrosynthesis_*.png #   逆合成分解スキーム図
    ├── 01_synthesis_*.png      #   前向き合成スキーム図 (AiZynthFinder)
    ├── ...
    ├── retrosynthesis_detail.json  # 詳細解析データ
    └── retrosynthesis_report.html  # HTML レポート (ブラウザで開く)
```

---

## 最終選抜結果 (Result_Best)

### 選抜基準

| 基準 | 内容 |
|------|------|
| DrugLike | Lipinski Ro5 + Veber ルール充足 |
| スコアカットオフ | ≤ -5.0 kcal/mol |
| ランキング | ドッキングスコア上位 + LE 上位の統合 Top 5 |

### LE (Ligand Efficiency) の見方

```
LE = |ドッキングスコア| / 重原子数 (HAC)

◎ LE ≥ 0.4  : 優秀 (フラグメント〜低分子の理想域)
○ LE ≥ 0.3  : 良好 (低分子医薬品の目安)
△ LE ≥ 0.2  : 許容範囲
✕ LE < 0.2  : 非効率 (大きすぎ or 弱すぎ)
```

### PyMOL での可視化

```python
# 受容体 + 最良候補の可視化
load results/docking/receptor.pdb
load Result_Best/01_score*.sdf

# 受容体を半透明表面表示
show surface, receptor
set transparency, 0.5, receptor
# 候補分子をスティック表示
show sticks, 01_score*
```

### 注意事項

- smina スコアは**同一受容体での候補間比較**に有効です。異なるターゲット間のスコア比較はできません
- SA Score / QED は計算予測値です。実際の合成可能性は合成化学者の判断が必要です
- PAINS / BRENK に引っかかっても偽陽性の場合があります。個別確認を推奨します
- **実験的検証** (IC50, SPR, ITC 等) による確認が最終的に必要です

---

## ファイル構成

```
Peptide_to_small_molecule/
│
│  # メインパイプライン
├── pipeline.py                 # 全 STEP 統合パイプライン
├── collect_best.py             # 最終候補選抜 + 合成スキーム図生成
│
│  # STEP 1: 相互作用解析
├── analyze_interactions.py     # 距離ベース接触検出
├── analyze_sasa.py             # ΔSASA 解析
├── analyze_prodigy.py          # PRODIGY 結合親和性予測
│
│  # STEP 2: ファーマコフォア
├── extract_pharmacophore.py    # ファーマコフォア特徴点抽出
│
│  # STEP 3: 低分子設計
├── design_small_molecule.py    # ペプチドミメティクス設計
├── pharmacophore_bridge.py     # ファーマコフォアブリッジ設計
├── rigid_scaffold_design.py    # 剛直スキャフォールド設計
├── generate_pocket_molecules.py# ポケット相補分子生成
├── linker_library.py           # リンカーライブラリ管理
│
│  # STEP 3b/3c: 合成評価
├── synthesizability.py         # SA Score / QED / PAINS / BRENK
├── retrosynthesis.py           # 逆合成解析 (BRICS / RECAP / AiZynthFinder)
│
│  # STEP 4: ドッキング
├── dock_with_smina.py          # smina ドッキング + リスコアリング
├── compare_cyclic_peptide.py   # 環状ペプチド比較
├── dock_cyclic_adcp.py         # AutoDock CrankPep
├── dock_smol_adfr.py           # AutoDockFR
│
│  # ユーティリティ
├── visualize.py                # 可視化
├── plot_utils.py               # matplotlib 設定
├── utils/                      # 共有モジュール
│   ├── drug_likeness.py        #   Ro5 + Veber ルール
│   ├── ligand_efficiency.py    #   LE 計算
│   └── residue_defs.py         #   アミノ酸定義・SMILES
│
│  # データ
├── Protein_Peptide.pdb         # サンプル入力 (タンパク質-ペプチド複合体)
├── linker_db.json              # リンカー DB (初回実行時に自動生成)
├── Linkers_from_FEgrow/        # FEgrow リンカーライブラリ
│   ├── library.sdf
│   └── smiles.txt
│
│  # セットアップ
├── adcpsuite_micromamba.sh     # ADCP/ADFR 環境セットアップ
├── LICENSE                     # MIT License
└── README.md                   # このファイル
```

---

## トラブルシューティング

### Q: `ModuleNotFoundError: No module named 'rdkit'`

conda 環境が有効化されていません:
```bash
conda activate peptide_pipeline
```

### Q: `smina: command not found`

smina がインストールされていないか、パスが通っていません:
```bash
# conda-forge からインストール
conda install -c conda-forge smina -y

# または静的バイナリをダウンロード
wget https://sourceforge.net/projects/smina/files/smina.static/download -O smina.static
chmod +x smina.static
export PATH="$PWD:$PATH"
```

### Q: `RuntimeError: Blosc2 library not found` (AiZynthFinder 使用時)

pytables を conda-forge から再インストールしてください:
```bash
conda install -c conda-forge pytables -y
```

### Q: AiZynthFinder のモデルデータが見つからない

`download_public_data` でモデルをダウンロードしてください:
```bash
download_public_data aizynthfinder_data
```

デフォルトでは以下の場所を自動検索します:
1. `./aizynthfinder_data/config.yml`
2. `~/aizynthfinder_data/config.yml`

別の場所に置いた場合は `--aizynthfinder-config` で指定してください。

### Q: ドッキングが遅い

`--exhaustiveness` を下げるとドッキングが速くなります (精度は下がります):
```bash
python pipeline.py input.pdb --exhaustiveness 4
```

### Q: 自分の PDB でエラーが出る

- PDB にタンパク質とペプチドが**別チェーン**で含まれているか確認してください
- 水分子やリガンドが含まれていても自動で除外されますが、問題がある場合は事前に除去してください
- チェーン ID を `--protein-chain` と `--peptide-chain` で明示的に指定してみてください

---

## ライセンス

MIT License - 詳細は [LICENSE](LICENSE) をご覧ください。
