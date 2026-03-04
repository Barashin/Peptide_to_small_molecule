"""
plot_utils.py
=============
matplotlib の日本語フォント設定ユーティリティ。
各描画モジュールの先頭でインポートして使用する。
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 日本語フォントを優先順位順に設定 (Linux / macOS 両対応)
plt.rcParams["font.family"] = ["Noto Sans CJK JP",
                                "Hiragino Sans", "Hiragino Maru Gothic Pro",
                                "Arial Unicode MS", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
