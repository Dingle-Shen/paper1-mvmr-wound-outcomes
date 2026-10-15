# -*- coding: utf-8 -*-
"""生成修正后的单暴露对照文件，供图 1 使用。

论文一的图 1 使用 mvmr_docA_r9local_summary.csv 的 uni_beta 列。
该列源自错误的加权公式，需替换为标准 IVW 的结果。

此处从 mediation_ivs_5outcomes.tsv（352 IV × 5 结局的暴露侧与结局侧效应，
随仓库发布）按修正公式重算单暴露估计，写回一个新的 summary 副本，
并把图表脚本的输出目录指定回 figures/。
"""
import os
import csv
import math
import numpy as np
import pandas as pd

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
RES = os.path.join(BASE, "results")
SRC = os.path.join(RES, "mvmr_docA_r9local_summary.csv")
IVS = os.path.join(RES, "mediation_ivs_5outcomes.tsv")
OUT = os.path.join(RES, "mvmr_docA_r9local_summary_corrected.csv")

d = pd.read_csv(IVS, sep="\t")
rows = []
for _, r in d.groupby("outcome"):
    pass

groups = {}
for oc, s in d.groupby("outcome"):
    groups[oc] = s

summary = pd.read_csv(SRC)
corrected = summary.copy()

n_fixed = 0
n_skip = 0
for idx, row in corrected.iterrows():
    oc = row["outcome"]
    exp = row["exposure"]
    s = groups.get(oc)
    if s is None or len(s) == 0:
        n_skip += 1
        corrected.at[idx, "uni_beta"] = np.nan
        continue
    col = "beta_" + exp
    if col not in s.columns:
        n_skip += 1
        corrected.at[idx, "uni_beta"] = np.nan
        continue
    B = s[col].values
    y = s["beta_out"].values
    sy = s["se_out"].values
    w = 1.0 / sy ** 2
    wj = w * B ** 2
    if np.sum(wj) <= 0:
        n_skip += 1
        continue
    ratio = float(np.sum(wj * (y / B)) / np.sum(wj))
    corrected.at[idx, "uni_beta"] = round(ratio, 4)
    n_fixed += 1

corrected.to_csv(OUT, index=False, encoding="utf-8-sig")
print("修正 %d 行，跳过 %d 行（无对应缓存）" % (n_fixed, n_skip))
print("已保存:", OUT)

# 打印对照
print("\n=== 对照：原值 vs 修正值 ===")
for idx, row in corrected.iterrows():
    oc = row["outcome"]
    exp = row["exposure"]
    old = summary.at[idx, "uni_beta"]
    new = row["uni_beta"]
    if new is None or (isinstance(new, float) and np.isnan(new)):
        print("  %-8s %-6s 原 %+10s -> 修正（缓存缺失，待重跑 R9）" % (oc, exp, old))
    else:
        print("  %-8s %-6s 原 %+10s -> 修正 %+10.4f" % (oc, exp, old, float(new)))
