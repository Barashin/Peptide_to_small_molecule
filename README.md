# Peptide to Small Molecule

ペプチド-タンパク質複合体の結合を保ちながら、ペプチドを低分子に変換するツールです。

---

## 1. 環境構築

```bash
# conda 環境を作成
conda create -n peptide_pipeline python=3.11 rdkit -c conda-forge -y
conda activate peptide_pipeline

# 必要パッケージをインストール
pip install biopython numpy scipy matplotlib pillow

# ドッキングエンジン (smina) をインストール
conda install -c conda-forge smina -y
```

### オプション: AI 逆合成 (AiZynthFinder)

より正確な合成ルートを提案させたい場合:

```bash
pip install aizynthfinder[all]
conda install -c conda-forge pytables -y
download_public_data aizynthfinder_data
```

---

## 2. 必要な入力

**PDB ファイル 1 つだけ** — タンパク質とペプチドが別チェーンで入った複合体構造。

サンプルファイル `Protein_Peptide.pdb` が同梱されているので、そのまま試せます。

---

## 3. 実行

```bash
conda activate peptide_pipeline

# サンプルで実行
python pipeline.py Protein_Peptide.pdb

# 自分の PDB で実行 (チェーン指定)
python pipeline.py my_complex.pdb --protein-chain A --peptide-chain B

# AI 逆合成を使う場合
python pipeline.py input.pdb --use-aizynthfinder
```

パイプライン完了後、上位候補の選抜と合成スキーム図を生成:

```bash
python collect_best.py
```

---

## 4. 結果

### `results/` — 全解析結果

- `contacts.json` — 相互作用データ
- `residue_scores.png` — 残基重要度グラフ
- `candidate_ligands.sdf` — 候補分子
- `docking/` — ドッキング結果 (スコア一覧 CSV + ポーズ SDF)

### `Result_Best/` — 上位 5 候補

- `*.sdf` — 選抜分子のドッキングポーズ (PyMOL で可視化)
- `summary.csv` — スコア・物性データ一覧
- `retrosynthesis/` — 合成スキーム図 (PNG) + HTML レポート

---

## ライセンス

MIT License
