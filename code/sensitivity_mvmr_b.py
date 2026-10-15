# -*- coding: utf-8 -*-
"""敏感性分析（修正版）：b_mvmr = 352 IV 全集 3 暴露 MVMR 的 CRP 直接效应（docA 同口径）。
对比单暴露 b 与 MVMR b 的中介比例。c 用 docA uni_beta（BMI 总效应权威口径）。
输出 sensitivity_mediation.tsv
"""
import os, sys, csv, math
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
import scipy.stats as st

RESULTS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
EXPOSURES = ["T2D", "BMI", "CRP"]
OUTCOMES = ["下肢溃疡", "皮肤感染", "蜂窝织炎", "慢性皮肤溃疡", "压疮"]


def load_raw():
    raw = {}
    with open(os.path.join(RESULTS, "mediation_exposure_raw.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d = {}
            for e in EXPOSURES:
                try:
                    d[e] = {"beta": float(r[f"beta_{e}"]), "se": float(r[f"se_{e}"]),
                            "ea": r[f"ea_{e}"], "nea": r[f"nea_{e}"]}
                except (TypeError, ValueError):
                    d[e] = None
            raw[r["rsid"]] = d
    return raw


def load_out():
    out = {}
    with open(os.path.join(RESULTS, "mediation_ivs_5outcomes.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[(r["outcome"], r["rsid"])] = {"beta": float(r["beta_out"]), "se": float(r["se_out"])}
    return out


def harmonize_to(b, se, ea, nea, ref_e, ref_a):
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


def mvmr_ivw(B, y, w):
    W = np.diag(w)
    XtWX = B.T @ W @ B
    XtWX_inv = np.linalg.pinv(XtWX)
    theta = XtWX_inv @ (B.T @ W @ y)
    se = np.sqrt(np.diag(XtWX_inv))
    z = theta / se
    p = 2 * (1 - st.norm.cdf(np.abs(z)))
    return theta, se, p


def main():
    raw = load_raw()
    outd = load_out()
    import gzip
    FULL = os.environ.get("MR_FINNGEN_R9_DIR", os.path.join(RESULTS, "finngen_full"))
    FNAMES = {
        "下肢溃疡": "finngen_R9_L12_ULCERLOWLIMB.gz",
        "皮肤感染": "finngen_R9_L12_INFECT_SKIN.gz",
        "蜂窝织炎": "finngen_R9_L12_CELLULITIS.gz",
        "压疮": "finngen_R9_L12_DECUBITANSULCERANDPRESSURE.sumstats.txt.gz",
        "慢性皮肤溃疡": "finngen_R9_L12_CHRONICULCEROFSKIN.gz",
    }
    want = set(raw.keys())
    alleles = {}
    for cn, fn in FNAMES.items():
        rec = {}
        with gzip.open(os.path.join(FULL, fn), "rt", encoding="utf-8", errors="replace") as f:
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

    # docA 权威值：CRP MVMR 直接效应、BMI uni 总效应
    doc = {}
    with open(os.path.join(RESULTS, "mvmr_docA_r9local_summary_corrected.csv"), encoding="utf-8") as f:
        for r in csv.DictReader(f):
            doc[(r["\ufeffoutcome"], r["exposure"])] = r

    a = 0.3378
    a_se = 0.0087
    print(f"{'结局':<8}{'a':>7}{'b_single':>9}{'b_mvmr(docA)':>13}{'se':>7}"
          f"  |  prop_single:{'':<2}{'prop_mvmr':>9}  |  ab_single:{'':<2}{'ab_mvmr':>9}")
    rows = []
    for oc in OUTCOMES:
        oa_oc = alleles[oc]
        # 352 全集的 3 暴露 MVMR（以 CRP 为参考等位）
        B_list, y_list, w_list = [], [], []
        for rs, d in raw.items():
            rC, rB, rT = d.get("CRP"), d.get("BMI"), d.get("T2D")
            o = outd.get((oc, rs))
            if not rC or not rB or not rT or not o:
                continue
            oref, oalt = oa_oc.get(rs, (None, None))
            if not oref or not oalt:
                continue
            hO = harmonize_to(o["beta"], o["se"], oalt, oref, rC["ea"], rC["nea"])
            if hO is None:
                continue
            hB = harmonize_to(rB["beta"], rB["se"], rB["ea"], rB["nea"], rC["ea"], rC["nea"])
            if hB is None:
                continue
            hT = harmonize_to(rT["beta"], rT["se"], rT["ea"], rT["nea"], rC["ea"], rC["nea"])
            if hT is None:
                continue
            B_list.append([rC["beta"], hB[0], hT[0]])
            y_list.append(hO[0])
            w_list.append(1.0 / o["se"] ** 2)
        B = np.array(B_list); y = np.array(y_list); w = np.array(w_list)
        theta, se, p = mvmr_ivw(B, y, w)
        b_mvmr, se_mvmr = theta[0], se[0]  # CRP 直接效应
        # 单暴露 b（CRP 工具）
        b1_list, b1_se = [], []
        for rs, d in raw.items():
            rC = d.get("CRP")
            o = outd.get((oc, rs))
            if not rC or not o:
                continue
            if abs(rC["beta"]) / rC["se"] < math.sqrt(10):
                continue
            oref, oalt = oa_oc.get(rs, (None, None))
            if not oref or not oalt:
                continue
            hO = harmonize_to(o["beta"], o["se"], oalt, oref, rC["ea"], rC["nea"])
            if hO is None:
                continue
            b1_list.append((rC["beta"], hO[0]))
            b1_se.append(o["se"])
        ratio = np.array([t[1] / t[0] for t in b1_list if t[0] != 0])
        se_r = np.array(b1_se)[np.array([t[0] for t in b1_list]) != 0] / np.abs(np.array([t[0] for t in b1_list])[np.array([t[0] for t in b1_list]) != 0])
        ww = 1.0 / se_r ** 2
        b_single = float(np.sum(ww * ratio) / np.sum(ww))

        # docA c（BMI uni 总效应）
        c_val = float(doc[(oc, "BMI")]["uni_beta"])
        ab_s, ab_m = a * b_single, a * b_mvmr
        prop_s, prop_m = ab_s / c_val, ab_m / c_val
        rows.append({"outcome": oc, "a": round(a, 5), "b_single": round(b_single, 5),
                     "b_mvmr": round(b_mvmr, 5), "b_mvmr_se": round(se_mvmr, 5),
                     "ab_single": round(ab_s, 5), "ab_mvmr": round(ab_m, 5),
                     "prop_single": round(prop_s, 4), "prop_mvmr": round(prop_m, 4),
                     "c_docA": round(c_val, 5)})
        print(f"{oc:<8}{a:>7.4f}{b_single:>9.4f}{b_mvmr:>13.4f}{se_mvmr:>7.4f}"
              f"  |  {prop_s:>10.2%}{prop_m:>9.2%}  |  {ab_s:>10.4f}{ab_m:>9.4f}")

    with open(os.path.join(RESULTS, "sensitivity_mediation.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"写出 -> results/sensitivity_mediation.tsv")


if __name__ == "__main__":
    main()
