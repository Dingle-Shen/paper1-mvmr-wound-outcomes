# -*- coding: utf-8 -*-
"""论文一：多效性检验与单暴露审计。

产出三项：
  1) MVMR-Egger 截距与各暴露系数   -> results/mvmr_pleiotropy_egger.tsv
  2) 异质性 Q 与离群变异归因        -> results/mvmr_outlier_attribution.tsv
  3) 单暴露估计量审计（修正前后）    -> results/mvmr_singlevar_audit.csv

输入：
  results/mvmr_localr9_ivs.tsv                  逐 SNP 暴露与结局效应（等位已协调）
  results/mvmr_docA_r9local_summary_corrected.csv  docA 主口径 MVMR 与单暴露结果

说明：
  - Egger 截距在同时纳入三暴露的加权回归中估计（权重 1/se_out^2）
  - Q 为 MVMR 残差加权的 Cochran 统计量，自由度 = 工具数 - 暴露数
  - 离群归因剔除对 Q 贡献最大的 5 个变异后重估
  - 单暴露审计对比三个口径：多变量估计、修正前加权、标准逆方差加权

用法：python run_pleiotropy_sensitivity.py
"""
import os
import sys
import csv

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scipy_shim  # noqa: F401
import numpy as np  # noqa: E402
import scipy.stats as st  # noqa: E402

RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))
EXPOSURES = ["T2D", "BMI", "CRP"]
N_TOP = 5


def load_ivs():
    """逐结局的逐 SNP 暴露与结局效应。

    用 mediation_ivs_5outcomes.tsv 而非 mvmr_localr9_ivs.tsv：后者缺慢性皮肤溃疡
    （仅 4 结局 × 352 = 1408 行），前者为 5 结局 × 352 = 1760 行。两者在共有的
    4 个结局上数值一致（暴露效应由同一来源复制）。
    """
    d = {}
    with open(os.path.join(RESULTS, "mediation_ivs_5outcomes.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d.setdefault(r["outcome"], []).append(r)
    return d


def load_docA():
    """读 docA 汇总。

    mvmr_beta 与 ivw 标准值取自 corrected 文件；
    uni_beta_现有 必须取自**修正前**的 mvmr_docA_r9local_summary.csv——
    该列源自错误加权公式，正是审计要对照的对象。若误取 corrected 文件，
    对照会退化为自比自，审计失效。
    """
    doc = {}
    for tag, fname in (("cur", "mvmr_docA_r9local_summary_corrected.csv"),
                       ("old", "mvmr_docA_r9local_summary.csv")):
        path = os.path.join(RESULTS, fname)
        if not os.path.exists(path):
            raise FileNotFoundError("缺少 docA 文件：%s" % path)
        with open(path, encoding="utf-8") as f:
            for r in csv.DictReader(f):
                key = r.get("outcome") or r.get("\ufeffoutcome")
                e = doc.setdefault((key, r["exposure"]), {})
                e["mvmr_beta"] = float(r["mvmr_beta"])
                e[tag + "_uni"] = float(r["uni_beta"])
    return doc


def wls(X, y, w):
    """加权最小二乘：返回 (系数, 标准误, P)"""
    W = np.diag(w)
    XtWX_inv = np.linalg.pinv(X.T @ W @ X)
    beta = XtWX_inv @ (X.T @ W @ y)
    se = np.sqrt(np.diag(XtWX_inv))
    p = 2 * (1 - st.norm.cdf(np.abs(beta / se)))
    return beta, se, p


def build(rows, keep=None):
    """构造 B(n×3)、y、w"""
    B, y, w = [], [], []
    for i, r in enumerate(rows):
        if keep is not None and not keep[i]:
            continue
        B.append([float(r["beta_" + e]) for e in EXPOSURES])
        y.append(float(r["beta_out"]))
        w.append(1.0 / float(r["se_out"]) ** 2)
    return np.array(B), np.array(y), np.array(w)


def main():
    ivs = load_ivs()
    doc = load_docA()
    out = ivs["下肢溃疡"][0]["outcome"]  # 占位，仅取结局名集合
    outcomes = sorted(ivs.keys())

    egger_rows, out_rows, audit_rows = [], [], []
    for oc in outcomes:
        rows = ivs[oc]
        B, y, w = build(rows)

        # 1) MVMR 主效应（无截距）
        theta, se_t, p_t = wls(B, y, w)
        resid = y - B @ theta
        Q_full = float(np.sum(w * resid ** 2))
        df_full = len(rows) - len(EXPOSURES)
        Qp_full = float(st.chi2.sf(Q_full, df_full))

        # 2) MVMR-Egger：加截距列
        X1 = np.hstack([B, np.ones((len(rows), 1))])
        beta_e, se_e, p_e = wls(X1, y, w)
        egger_rows.append({
            "outcome": oc, "n_iv": len(rows),
            "egger_intercept": round(float(beta_e[3]), 5),
            "egger_se": round(float(se_e[3]), 5),
            "egger_p": round(float(p_e[3]), 3),
            "egger_T2D": round(float(beta_e[0]), 4),
            "egger_BMI": round(float(beta_e[1]), 4),
            "egger_CRP": round(float(beta_e[2]), 4),
        })

        # 3) 离群归因：按对 Q 的贡献取前 5
        contrib = w * resid ** 2
        order = np.argsort(-contrib)
        top_idx = set(order[:N_TOP].tolist())
        top_snps = ";".join(rows[i]["rsid"] for i in order[:N_TOP])
        keep = [i not in top_idx for i in range(len(rows))]
        Be, ye, we = build(rows, keep)
        theta_e, _, _ = wls(Be, ye, we)
        resid_e = ye - Be @ theta_e
        Q_excl = float(np.sum(we * resid_e ** 2))
        df_e = len(ye) - len(EXPOSURES)
        Qp_excl = float(st.chi2.sf(Q_excl, df_e))
        out_rows.append({
            "outcome": oc, "n_iv": len(rows),
            "Q_full": round(Q_full, 1), "Q_p_full": round(Qp_full, 4),
            "Q_excl5": round(Q_excl, 1), "Q_p_excl5": round(Qp_excl, 3),
            "top5_snps": top_snps,
            "excl5_T2D": round(float(theta_e[0]), 4),
            "excl5_BMI": round(float(theta_e[1]), 4),
            "excl5_CRP": round(float(theta_e[2]), 4),
            "full_T2D": round(float(theta[0]), 4),
            "full_BMI": round(float(theta[1]), 4),
            "full_CRP": round(float(theta[2]), 4),
        })

        # 4) 单暴露审计
        for k, e in enumerate(EXPOSURES):
            Bc = B[:, k]
            m = Bc != 0
            ratio = y[m] / Bc[m]
            wj = w[m] * Bc[m] ** 2          # 标准 IVW 权重
            ivw = float(np.sum(wj * ratio) / np.sum(wj))
            r = doc.get((oc, e), {})
            mv = r.get("mvmr_beta", float("nan"))
            old = r.get("old_uni", float("nan"))     # 修正前（错误加权）
            new = r.get("cur_uni", float("nan"))     # 修正后（标准 IVW，docA 权威列）
            audit_rows.append({
                "outcome": oc, "exposure": e, "n_iv": len(rows),
                "mvmr_beta": mv,
                "uni_beta_现有": old,
                "ivw_beta_标准": round(ivw, 4),
                "方向_mvmr": "+" if mv > 0 else "-",
                "方向_现有": "+" if old > 0 else "-",
                "方向_标准IVW": "+" if ivw > 0 else "-",
                "翻转_现有": "是" if (mv > 0) != (old > 0) else "否",
                "翻转_标准": "是" if (mv > 0) != (ivw > 0) else "否",
                "docA_修正后": new,
                "重算_与docA修正后差": round(abs(ivw - new), 5) if new == new else "",
            })

    # 统一写出（全部制表符分隔，修正原 .tsv 用了逗号的问题）
    def dump(name, rows, fields):
        p = os.path.join(RESULTS, name)
        with open(p, "w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader(); w.writerows(rows)
        print("写出 -> %s" % p)

    dump("mvmr_pleiotropy_egger.tsv", egger_rows,
         ["outcome", "n_iv", "egger_intercept", "egger_se", "egger_p",
          "egger_T2D", "egger_BMI", "egger_CRP"])
    dump("mvmr_outlier_attribution.tsv", out_rows,
         ["outcome", "n_iv", "Q_full", "Q_p_full", "Q_excl5", "Q_p_excl5", "top5_snps",
          "excl5_T2D", "excl5_BMI", "excl5_CRP", "full_T2D", "full_BMI", "full_CRP"])
    # 审计表保留 csv 扩展名（与正文引用一致），同为逗号分隔
    p = os.path.join(RESULTS, "mvmr_singlevar_audit.csv")
    with open(p, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(audit_rows[0].keys()))
        w.writeheader(); w.writerows(audit_rows)
    print("写出 -> %s" % p)

    print()
    print("=== Egger 截距 ===")
    for r in egger_rows:
        print("  %-8s intercept=%+.5f (SE %.5f, P=%.3f)" %
              (r["outcome"], r["egger_intercept"], r["egger_se"], r["egger_p"]))
    print()
    print("=== Q 与剔除后 ===")
    for r in out_rows:
        print("  %-8s Q %.1f (P=%.4f) -> 剔除5后 %.1f (P=%.3f)" %
              (r["outcome"], r["Q_full"], r["Q_p_full"], r["Q_excl5"], r["Q_p_excl5"]))
    print()
    nflip = sum(1 for r in audit_rows if r["翻转_现有"] == "是")
    nflip_s = sum(1 for r in audit_rows if r["翻转_标准"] == "是")
    print("单暴露审计：修正前方向背离 %d / 15，标准 IVW 方向背离 %d / 15" %
          (nflip, nflip_s))


if __name__ == "__main__":
    main()
