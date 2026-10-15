# -*- coding: utf-8 -*-
"""论文一：统计功效与可检出最小效应量（正文表 3）。

在 α=0.05（双侧）、功效 80% 的条件下，可检出的最小效应量为观测标准误的
z 倍，其中 z = z_{1-α/2} + z_{1-β} = 1.959964 + 0.841621 = 2.801585。
结局为二分类，故以 OR 表示：OR_MDE = exp(z * SE)。

SE 取多变量估计的标准误（与正文表 2 同口径）。

输入：results/mvmr_docA_r9local_summary_corrected.csv
输出：results/mde_power.tsv
      （outcome, n_cases, exposure, beta, se, mde_beta, mde_or）
      并打印正文表 3 的 OR 矩阵

用法：python run_mde_power.py
"""
import io
import csv
import math
import os
import sys

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))

# z_{1-0.05/2} + z_{1-0.20}
Z = 1.959963985 + 0.841621234

EXPOSURES = ["T2D", "BMI", "CRP"]
EXP_LABEL = {"T2D": "类型 2 糖尿病", "BMI": "体重指数", "CRP": "C 反应蛋白"}
OUTCOMES = ["皮肤感染", "蜂窝织炎", "下肢溃疡", "慢性皮肤溃疡", "压疮"]
N_CASES = {"皮肤感染": 21308, "蜂窝织炎": 6683, "下肢溃疡": 3784,
           "慢性皮肤溃疡": 1840, "压疮": 1479}


def main():
    src = os.path.join(RESULTS, "mvmr_docA_r9local_summary_corrected.csv")
    with open(src, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise SystemExit("输入为空：%s" % src)
    key0 = list(rows[0].keys())[0].replace("\ufeff", "")
    for r in rows:
        r[key0] = r.get(key0) or r.get("\ufeff" + key0)

    out = []
    grid = {}
    for oc in OUTCOMES:
        for e in EXPOSURES:
            m = [r for r in rows if (r.get("outcome") or r.get(key0)) == oc
                 and r["exposure"] == e]
            if not m:
                continue
            r = m[0]
            beta = float(r["mvmr_beta"])
            se = float(r["mvmr_se"])
            mde_beta = Z * se
            mde_or = math.exp(mde_beta)
            grid[(oc, e)] = mde_or
            out.append({
                "outcome": oc, "n_cases": N_CASES[oc], "exposure": e,
                "beta": beta, "se": se,
                "mde_beta": round(mde_beta, 4), "mde_or": round(mde_or, 3),
            })

    dst = os.path.join(RESULTS, "mde_power.tsv")
    with open(dst, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(out)
    print("z = %.6f" % Z)
    print("写出 -> %s" % dst)
    print()

    # 正文表 3 的呈现：行＝结局，列＝暴露，值＝OR
    print("表 3  各结局可检出的最小效应量（OR，α=0.05，功效 80%）")
    print("%-16s %9s %10s %10s %10s" % ("结局", "病例数", "糖尿病", "体重指数", "CRP"))
    for oc in OUTCOMES:
        vals = ["%.2f" % grid[(oc, e)] if (oc, e) in grid else "—" for e in EXPOSURES]
        print("%-16s %9d %10s %10s %10s" % (oc, N_CASES[oc], *vals))
    print()
    print("读法：低于该值的效应不能被本设计排除，只能表述为未检出。")


if __name__ == "__main__":
    main()
