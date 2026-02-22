# Result_Best — 選抜分子一覧


全ドッキング結果から **スコア上位 + Ligand Efficiency 上位** を選抜した最終候補。


## 選抜基準
| 基準 | 内容 |
|------|------|
| DrugLike | Lipinski Ro5 + Veber 充足 |
| スコアカットオフ | ≤ -5.0 kcal/mol |
| 統合上位 | Top 5 (スコア + LE の統合ランキング) |

## コンフォメーション適合性指標
| 指標 | 範囲 | 意味 |
|------|------|------|
| Conformance Rate | 0〜1 | 低エネルギーコンフォマーのうちファーマコフォア距離が許容範囲内の割合 (高い = 良好) |
| Rotatable Between | 0〜 | ファーマコフォア間の回転可能結合数 (少ない = 剛直) |

## 合成容易性指標
| 指標 | 範囲 | 意味 |
|------|------|------|
| SA Score | 1〜10 | 合成容易性 (低い = 合成しやすい, ≤6 推奨) |
| QED | 0〜1 | 薬剤らしさ (高い = 良好, ≥0.5 推奨) |
| PAINS | OK/NG | アッセイ干渉構造なし = OK |
| BRENK | OK/NG | 構造アラートなし = OK |

## LE (Ligand Efficiency) の見方
```
LE = |smina Affinity| / 重原子数 (HAC)

◎ LE ≥ 0.4  : 優秀 (フラグメント〜低分子の理想域)
○ LE ≥ 0.3  : 良好 (低分子医薬品の目安)
△ LE ≥ 0.2  : 許容範囲
✕ LE < 0.2  : 非効率 (大きすぎ or 弱すぎ)
```

## 選抜分子一覧

| 順位 | スコア (kcal/mol) | HAC | LE | LE評価 | SA | QED | P | B | Conf | Rot | 分子名 |
|------|-----------------|-----|-----|--------|-----|------|---|---|------|-----|--------|
| 1 | -7.400 | 26 | 0.2846 | △ 許容 | 4.0 | 0.450 | OK | NG | - | - | reduced_peptide_TRP_ALA_PRO_direct |
| 2 | -6.100 | 16 | 0.3812 | ○ 良好 | 3.9 | 0.781 | OK | OK | - | - | direct_link_TRP_ALA_PRO |
| 3 | -6.100 | 22 | 0.2773 | △ 許容 | 5.6 | 0.782 | OK | OK | - | - | benzene_tri_TRP_ALA_PRO |
| 4 | -6.000 | 11 | 0.5455 | ◎ 優秀 | 1.9 | 0.635 | OK | OK | - | - | direct_link_TRP_ALA_direct |
| 5 | -5.800 | 34 | 0.1706 | ✕ 非効率 | 3.9 | 0.203 | OK | NG | - | - | reduced_peptide_TRP_ALA_PRO_C4 |

## ファイル構成
```
Result_Best/
├── 01_score*_<name>.sdf        # ドッキングポーズ SDF (PyMOL/VMD で可視化可)
├── ...                         # 02〜N 番まで同様
├── retrosynthesis/             # 逆合成解析結果
│   ├── 01_retrosynthesis_*.png # 合成スキーム図
│   ├── ...                     # 各分子ごとの図
│   └── retrosynthesis_detail.json  # 全詳細データ (JSON)
├── summary.csv                 # 全選抜分子の詳細データ
└── README.md                   # 本ファイル
```

## PyMOL での可視化
```python
load results/docking/receptor.pdb
load Result_Best/01_score*.sdf
# 受容体を表面表示
show surface, receptor
set transparency, 0.5, receptor
```

## 注意事項
- smina (AutoDock Vina) スコアは**同一受容体での比較**に有効
- ペプチド vs 低分子の直接比較は LE を参照すること
- 実験的検証 (IC50, SPR, ITC 等) で確認が必要
