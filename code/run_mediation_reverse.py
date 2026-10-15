# -*- coding: utf-8 -*-
"""反向中介对照（方向敏感性分析）：CRP→BMI→结局。

目的：对中介方向设定（BMI→CRP→结局）作方向敏感性检验，非方向确认。
设计与解读规则见 reports/03_中介分析/70_中介分析设计定稿_四个决策点与整合判断.md
（§7.3 解读规则、§7.6 先诊断后对照）。

口径约束（与正向严格对称，见 70 §7.6）：
  - a、c 由暴露自身工具集估计（|Z|>=sqrt(10)）
  - b 主口径为 MVMR 联合工具集调整估计，敏感性口径为中介自身工具集
  - b 与 b_rev 同为 MVMR 调整口径；c 与 c_rev 同为单暴露总效应口径

输入：
  results/mediation_exposure_raw.tsv   352 IV × 三暴露原始 beta/se/ea/nea
  results/mediation_ivs_5outcomes.tsv  352 IV × 5 结局 beta_out/se_out
  results/sensitivity_mediation.tsv    正向结果（用于三分量对照）
  results/mediation_results_5outcomes.tsv 正向逐步结果（用于 a_fwd / c_fwd 核对）
  FinnGen R9 全量汇总统计（等位协调用；路径可用环境变量 MR_FINNGEN_R9_DIR 覆盖）

输出：
  results/mediation_reverse.tsv            反向中介逐结局结果
  results/mediation_reverse_diagnosis.tsv  三分量分解与档位判定（70 §7.6 三步）

用法：
  python run_mediation_reverse.py
"""
import os
import sys
import csv
import math
import gzip
import time

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import scipy_shim  # noqa: F401  注册最小 scipy 兼容层（本环境无 scipy）
import numpy as np  # noqa: E402
import scipy.stats as st  # noqa: E402

RESULTS = os.path.normpath(os.path.join(HERE, "..", "results"))
FULL = os.environ.get("MR_FINNGEN_R9_DIR", os.path.join(RESULTS, "finngen_full"))
CACHE = os.path.join(RESULTS, "_outcome_alleles_cache.tsv")

COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
EXPOSURES = ["T2D", "BMI", "CRP"]
OUTCOMES = ["下肢溃疡", "皮肤感染", "蜂窝织炎", "慢性皮肤溃疡", "压疮"]
FNAMES = {
    "下肢溃疡": "finngen_R9_L12_ULCERLOWLIMB.gz",
    "皮肤感染": "finngen_R9_L12_INFECT_SKIN.gz",
    "蜂窝织炎": "finngen_R9_L12_CELLULITIS.gz",
    "压疮": "finngen_R9_L12_DECUBITANSULCERANDPRESSURE.sumstats.txt.gz",
    "慢性皮肤溃疡": "finngen_R9_L12_CHRONICULCEROFSKIN.gz",
}
Z_MIN = math.sqrt(10.0)  # |Z|>=sqrt(10) 近似 F>=10，与正向一致


# ---------------------------------------------------------------- 数据载入
def load_raw():
    """{rsid: {exp: {beta, se, ea, nea}}}"""
    raw = {}
    with open(os.path.join(RESULTS, "mediation_exposure_raw.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d = {}
            for e in EXPOSURES:
                try:
                    d[e] = {"beta": float(r[f"beta_{e}"]), "se": float(r[f"se_{e}"]),
                            "ea": r[f"ea_{e}"], "nea": r[f"nea_{e}"]}
                except (TypeError, ValueError, KeyError):
                    d[e] = None
            raw[r["rsid"]] = d
    return raw


def load_out():
    """{(outcome, rsid): {beta, se}}"""
    out = {}
    with open(os.path.join(RESULTS, "mediation_ivs_5outcomes.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[(r["outcome"], r["rsid"])] = {"beta": float(r["beta_out"]),
                                              "se": float(r["se_out"])}
    return out


def load_outcome_alleles(want):
    """{outcome: {rsid: (ref, alt)}}。带本地缓存，避免每次重扫 3.8 GB。"""
    if os.path.exists(CACHE):
        cache = {}
        with open(CACHE, encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                cache.setdefault(r["outcome"], {})[r["rsid"]] = (r["ref"], r["alt"])
        if all(len(cache.get(o, {})) > 0 for o in OUTCOMES):
            print("  等位表命中缓存 -> %s" % os.path.basename(CACHE))
            return cache
    alleles, rows = {}, []
    for cn, fn in FNAMES.items():
        path = os.path.join(FULL, fn)
        if not os.path.exists(path):
            raise FileNotFoundError("缺少 FinnGen 文件：%s" % path)
        t0 = time.time()
        rec = {}
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                p = line.rstrip("\n").split("\t")
                if len(p) < 12 or not p[4]:
                    continue
                for r in p[4].split(";"):
                    r = r.strip()
                    if r in want:
                        rec[r] = (p[2], p[3])
                        break
        alleles[cn] = rec
        for rs, (rf, al) in rec.items():
            rows.append({"outcome": cn, "rsid": rs, "ref": rf, "alt": al})
        print("  结局等位 %-8s %3d/%d  %.0fs" % (cn, len(rec), len(want), time.time() - t0))
    with open(CACHE, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["outcome", "rsid", "ref", "alt"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    return alleles


# ---------------------------------------------------------------- 统计工具
def harmonize_to(b, se, ea, nea, ref_e, ref_a):
    """把效应协调到参考等位方向；等位冲突返回 None。"""
    if ea == ref_e:
        return b, se
    if ea == ref_a:
        return -b, se
    c = COMP.get(ea)
    if c == ref_e:
        return b, se
    if c == ref_a:
        return -b, se
    return None


def ivw(bx, by, sey):
    """Wald ratio + 逆方差加权（与正向脚本同实现）。"""
    bx = np.asarray(bx, float); by = np.asarray(by, float); sey = np.asarray(sey, float)
    m = ~(np.isnan(bx) | np.isnan(by) | np.isnan(sey) | (bx == 0))
    bx, by, sey = bx[m], by[m], sey[m]
    if len(bx) < 3:
        return None
    ratio = by / bx
    se_ratio = sey / np.abs(bx)
    w = 1.0 / se_ratio ** 2
    b_ivw = float(np.sum(w * ratio) / np.sum(w))
    se_ivw = float(np.sqrt(1.0 / np.sum(w)))
    q = float(np.sum(w * (ratio - b_ivw) ** 2))
    df = len(bx) - 1
    p_q = float(st.chi2.sf(q, df)) if df > 0 else 1.0
    z = b_ivw / se_ivw
    p = 2 * (1 - st.norm.cdf(abs(z)))
    return {"beta": b_ivw, "se": se_ivw, "p": p, "n_snp": len(bx), "Q": q, "Q_p": p_q}


def mvmr_ivw(B, y, w):
    """多变量逆方差加权（与 sensitivity_mvmr_b.py 同实现）。"""
    W = np.diag(w)
    XtWX_inv = np.linalg.pinv(B.T @ W @ B)
    theta = XtWX_inv @ (B.T @ W @ y)
    se = np.sqrt(np.diag(XtWX_inv))
    z = theta / se
    p = 2 * (1 - st.norm.cdf(np.abs(z)))
    return theta, se, p


def load_doca_uni():
    """docA 权威单暴露总效应：全 352 IV、权重 wj = (1/se_out^2)*B^2。

    正向主口径的 c 即取自此列（sensitivity_mvmr_b.py 的 c_docA）。
    为满足 70 §7.6 的口径对称要求，反向的 c_rev 必须同用此口径。
    返回 {(outcome, exposure): uni_beta}
    """
    doc = {}
    with open(os.path.join(RESULTS, "mvmr_docA_r9local_summary_corrected.csv"),
              encoding="utf-8") as f:
        for r in csv.DictReader(f):
            key = r.get("outcome") or r.get("\ufeffoutcome")
            doc[(key, r["exposure"])] = float(r["uni_beta"])
    return doc


def load_forward():
    """正向结果：sensitivity_mediation.tsv + mediation_results_5outcomes.tsv"""
    fwd = {}
    with open(os.path.join(RESULTS, "sensitivity_mediation.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            fwd[r["outcome"]] = {"a_fwd": float(r["a"]),
                                 "b_fwd_single": float(r["b_single"]),
                                 "b_fwd_mvmr": float(r["b_mvmr"]),
                                 "c_fwd": float(r["c_docA"]),
                                 "prop_fwd_single": float(r["prop_single"]),
                                 "prop_fwd_mvmr": float(r["prop_mvmr"])}
    # 用原始中介表核对 a 与 c（口径一致性自检）
    with open(os.path.join(RESULTS, "mediation_results_5outcomes.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            o = r["outcome"]
            if o in fwd:
                fwd[o]["a_raw"] = float(r["a_bmi_to_crp"])
                fwd[o]["c_raw"] = float(r["c_total"])
    return fwd


# ---------------------------------------------------------------- 主流程
def main():
    raw = load_raw()
    outd = load_out()
    fwd = load_forward()
    doca = load_doca_uni()
    print("IV: %d | 结局记录: %d | 正向结局: %d | docA uni: %d"
          % (len(raw), len(outd), len(fwd), len(doca)))
    alleles = load_outcome_alleles(set(raw.keys()))

    rows, diag = [], []
    for oc in OUTCOMES:
        oa = alleles[oc]

        # ---- a_rev: CRP → BMI（CRP 工具） ----
        bx_a, by_a, se_a = [], [], []
        for rs, d in raw.items():
            rC, rB = d.get("CRP"), d.get("BMI")
            if not rC or not rB:
                continue
            if abs(rC["beta"]) / rC["se"] < Z_MIN:
                continue
            h = harmonize_to(rB["beta"], rB["se"], rB["ea"], rB["nea"], rC["ea"], rC["nea"])
            if h is None:
                continue
            bx_a.append(rC["beta"]); by_a.append(h[0]); se_a.append(rB["se"])
        a_rev = ivw(bx_a, by_a, se_a)

        # ---- b_rev（敏感性口径，单暴露）: BMI → 结局（BMI 工具） ----
        bx_b, by_b, se_b = [], [], []
        for rs, d in raw.items():
            rB = d.get("BMI")
            o = outd.get((oc, rs))
            if not rB or not o:
                continue
            if abs(rB["beta"]) / rB["se"] < Z_MIN:
                continue
            oref, oalt = oa.get(rs, (None, None))
            if not oref or not oalt:
                continue
            h = harmonize_to(o["beta"], o["se"], oalt, oref, rB["ea"], rB["nea"])
            if h is None:
                continue
            bx_b.append(rB["beta"]); by_b.append(h[0]); se_b.append(o["se"])
        b_rev_single = ivw(bx_b, by_b, se_b)

        # ---- c_rev: CRP → 结局（CRP 工具，单暴露总效应口径） ----
        bx_c, by_c, se_c = [], [], []
        for rs, d in raw.items():
            rC = d.get("CRP")
            o = outd.get((oc, rs))
            if not rC or not o:
                continue
            if abs(rC["beta"]) / rC["se"] < Z_MIN:
                continue
            oref, oalt = oa.get(rs, (None, None))
            if not oref or not oalt:
                continue
            h = harmonize_to(o["beta"], o["se"], oalt, oref, rC["ea"], rC["nea"])
            if h is None:
                continue
            bx_c.append(rC["beta"]); by_c.append(h[0]); se_c.append(o["se"])
        c_rev = ivw(bx_c, by_c, se_c)

        # ---- b_rev（主口径，MVMR 联合工具集）: BMI 直接效应 ----
        # 参考等位取 CRP；MVMR 系数对等位参考方向不变
        B_list, y_list, w_list = [], [], []
        for rs, d in raw.items():
            rC, rB, rT = d.get("CRP"), d.get("BMI"), d.get("T2D")
            o = outd.get((oc, rs))
            if not rC or not rB or not rT or not o:
                continue
            oref, oalt = oa.get(rs, (None, None))
            if not oref or not oalt:
                continue
            hO = harmonize_to(o["beta"], o["se"], oalt, oref, rC["ea"], rC["nea"])
            hB = harmonize_to(rB["beta"], rB["se"], rB["ea"], rB["nea"], rC["ea"], rC["nea"])
            hT = harmonize_to(rT["beta"], rT["se"], rT["ea"], rT["nea"], rC["ea"], rC["nea"])
            if hO is None or hB is None or hT is None:
                continue
            B_list.append([rC["beta"], hB[0], hT[0]])
            y_list.append(hO[0])
            w_list.append(1.0 / o["se"] ** 2)
        theta, se_theta, p_theta = mvmr_ivw(np.array(B_list), np.array(y_list), np.array(w_list))
        b_rev_mvmr, se_b_rev_mvmr = float(theta[1]), float(se_theta[1])   # 索引 1 = BMI

        if not (a_rev and b_rev_single and c_rev):
            print("  %s: SNP 不足，跳过" % oc)
            continue

        a_r = a_rev["beta"]
        # c_rev 主口径 = docA 全 352 IV 单暴露总效应（与正向 c_docA 同口径）
        c_rev_doca = doca[(oc, "CRP")]
        c_rev_filt = c_rev["beta"]          # 敏感性：CRP 工具 |Z|>=sqrt(10) 子集
        ab_s = a_r * b_rev_single["beta"]
        se_ab_s = math.sqrt((a_r * b_rev_single["se"]) ** 2 +
                            (b_rev_single["beta"] * a_rev["se"]) ** 2)
        ab_m = a_r * b_rev_mvmr
        se_ab_m = math.sqrt((a_r * se_b_rev_mvmr) ** 2 +
                            (b_rev_mvmr * a_rev["se"]) ** 2)
        prop_s = ab_s / c_rev_doca if c_rev_doca else float("nan")
        prop_m = ab_m / c_rev_doca if c_rev_doca else float("nan")
        prop_s_filt = ab_s / c_rev_filt if c_rev_filt else float("nan")
        prop_m_filt = ab_m / c_rev_filt if c_rev_filt else float("nan")
        p_ab_s = 2 * (1 - st.norm.cdf(abs(ab_s / se_ab_s)))
        p_ab_m = 2 * (1 - st.norm.cdf(abs(ab_m / se_ab_m)))

        rows.append({
            "outcome": oc,
            "a_rev": round(a_r, 5), "a_rev_se": round(a_rev["se"], 5),
            "a_rev_p": "%.3g" % a_rev["p"], "a_rev_n": a_rev["n_snp"],
            "b_rev_single": round(b_rev_single["beta"], 5),
            "b_rev_single_se": round(b_rev_single["se"], 5),
            "b_rev_single_n": b_rev_single["n_snp"],
            "b_rev_mvmr": round(b_rev_mvmr, 5), "b_rev_mvmr_se": round(se_b_rev_mvmr, 5),
            "c_rev_docA": round(c_rev_doca, 5), "c_rev_docA_n": 352,
            "c_rev_filtered": round(c_rev_filt, 5), "c_rev_filtered_n": c_rev["n_snp"],
            "c_rev_p_filtered": "%.3g" % c_rev["p"],
            "ab_rev_single": round(ab_s, 5), "p_ab_rev_single": "%.3g" % p_ab_s,
            "ab_rev_mvmr": round(ab_m, 5), "p_ab_rev_mvmr": "%.3g" % p_ab_m,
            "prop_rev_single": round(prop_s, 4),
            "prop_rev_mvmr": round(prop_m, 4),
            "prop_rev_single_filt": round(prop_s_filt, 4),
            "prop_rev_mvmr_filt": round(prop_m_filt, 4),
        })

        # ---- 三分量诊断（70 §7.6 三步）----
        f = fwd[oc]
        r_a = a_r / f["a_fwd"]
        r_b = b_rev_mvmr / f["b_fwd_mvmr"]
        r_c = f["c_fwd"] / c_rev_doca
        ratio = prop_m / f["prop_fwd_mvmr"] if f["prop_fwd_mvmr"] else float("nan")
        # 档位判定（70 §7.3 修订表）：仅 a 项具判别力
        if ratio < 0.5:
            tier, verdict = "第一档", "方向设定得到数据支持"
        elif abs(ratio - 1) < 0.30:
            tier, verdict = "第二档", "数据无法区分方向"
        elif r_a <= 1:
            tier, verdict = "第三档第一行", "差异属效应强度不对称的数学产物，不构成方向证据"
        else:
            tier, verdict = "第三档第二行", "方向设定需重新审视"
        diag.append({
            "outcome": oc,
            "a_fwd": round(f["a_fwd"], 5), "a_rev": round(a_r, 5), "a_ratio": round(r_a, 4),
            "b_fwd_mvmr": round(f["b_fwd_mvmr"], 5), "b_rev_mvmr": round(b_rev_mvmr, 5),
            "b_ratio": round(r_b, 4),
            "c_fwd": round(f["c_fwd"], 5), "c_rev": round(c_rev_doca, 5),
            "c_ratio": round(r_c, 4),
            "ratio_product": round(r_a * r_b * r_c, 4),
            "prop_fwd_mvmr": round(f["prop_fwd_mvmr"], 4),
            "prop_rev_mvmr": round(prop_m, 4),
            "rev_over_fwd": round(ratio, 3),
            "tier": tier, "verdict": verdict,
        })

        print("  %-8s a_rev=%+.4f b_rev(MVMR)=%+.4f c_rev=%+.4f | prop: %.4f -> %.4f (%.2fx) [%s]"
              % (oc, a_r, b_rev_mvmr, c_rev_doca, f["prop_fwd_mvmr"], prop_m, ratio, tier))

    # ---- 写出 ----
    out1 = os.path.join(RESULTS, "mediation_reverse.tsv")
    with open(out1, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(rows)
    print("写出 -> results/mediation_reverse.tsv")

    out2 = os.path.join(RESULTS, "mediation_reverse_diagnosis.tsv")
    with open(out2, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(diag[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(diag)
    print("写出 -> results/mediation_reverse_diagnosis.tsv")

    # ---- 汇总判定 ----
    print()
    print("=== 三分量诊断汇总（a 项比值判别）===")
    for d in diag:
        print("  %-8s a比 %.4f | b比 %8.4f | c比 %6.4f | 乘积 %8.3f | %s | %s"
              % (d["outcome"], d["a_ratio"], d["b_ratio"], d["c_ratio"],
                 d["ratio_product"], d["tier"], d["verdict"]))
    all_a = all(d["a_ratio"] <= 1 for d in diag)
    print()
    print("全部结局 a 项比值 <= 1：%s" % ("是" if all_a else "否"))
    print("覆盖档位：%s" % sorted({d["tier"] for d in diag}))


if __name__ == "__main__":
    main()
