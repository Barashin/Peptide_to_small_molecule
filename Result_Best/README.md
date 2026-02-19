# Result_Best — 選抜分子一覧

詳細な選抜基準・分子一覧・LE評価は [メインREADME](../README.md#最終選抜結果-result_best) を参照してください。

## ファイル構成

```
Result_Best/
├── 01_score*_<name>.sdf   # ドッキングポーズ SDF (PyMOL/VMD で可視化可)
├── ...                    # 02〜N 番まで同様
├── summary.csv            # 全選抜分子の詳細データ
└── README.md              # 本ファイル
```

## 再生成

```bash
cd ..
python collect_best.py
```
