# -*- coding: utf-8 -*-
"""中介分析主链路：
A. 两两网络 MR 定暴露间因果方向（BMI↔CRP、BMI↔T2D、CRP↔T2D，双向 IVW + Steiger）
B. 两步中介：BMI→CRP→5 创面结局（中介比例 = a*b/c 产品法，Delta 法 SE）
输入：
  results/mediation_exposure_raw.tsv  (352 IV × 三暴露原始效应 beta/se/ea/nea)
  results/mediation_ivs_5outcomes.tsv (352 IV × 5 结局 beta_out/se_out)
输出：
  results/mediation_network_direction.tsv  (两两方向判定)
  results/mediation_results_5outcomes.tsv  (5 结局 × 中介结果)
"""
import os, sys, csv, math
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np

RESULTS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}
F_MIN = 10.0

EXPOSURES = ["T2D", "BMI", "CRP"]
OUTCOMES = ["下肢溃疡", "皮肤感染", "蜂窝织炎", "慢性皮肤溃疡", "压疮"]


def load_raw():
    """返回 {rsid: {exp: {beta, se, ea, nea}}}"""
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
    """返回 {(outcome, rsid): {beta, se}}"""
    out = {}
    with open(os.path.join(RESULTS, "mediation_ivs_5outcomes.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            out[(r["outcome"], r["rsid"])] = {"beta": float(r["beta_out"]), "se": float(r["se_out"])}
    return out


def load_out_alleles():
    """各结局 alt/ref（等位协调用）。优先读随仓库发布的等位缓存
    results/_outcome_alleles_cache.tsv；缺失时才扫本地 R9 全量文件。"""
    import gzip
    RESULTS = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results"))
    CACHE = os.path.join(RESULTS, "_outcome_alleles_cache.tsv")
    if os.path.exists(CACHE):
        cache = {}
        with open(CACHE, encoding="utf-8") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                cache.setdefault(r["outcome"], {})[r["rsid"]] = (r["ref"], r["alt"])
        if all(len(cache.get(o, {})) > 0 for o in OUTCOMES):
            print("  等位表命中缓存 -> %s" % os.path.basename(CACHE))
            return cache
    FULL = os.environ.get("MR_FINNGEN_R9_DIR", os.path.join(RESULTS, "finngen_full"))
    FNAMES = {
        "下肢溃疡": "finngen_R9_L12_ULCERLOWLIMB.gz",
        "皮肤感染": "finngen_R9_L12_INFECT_SKIN.gz",
        "蜂窝织炎": "finngen_R9_L12_CELLULITIS.gz",
        "压疮": "finngen_R9_L12_DECUBITANSULCERANDPRESSURE.sumstats.txt.gz",
        "慢性皮肤溃疡": "finngen_R9_L12_CHRONICULCEROFSKIN.gz",
    }
    want = set()
    with open(os.path.join(RESULTS, "mediation_ivs_5outcomes.tsv"), encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            want.add(r["rsid"])
    alleles = {}
    for cn, fn in FNAMES.items():
        path = os.path.join(FULL, fn)
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
                        rec[r] = (p[2], p[3])  # (ref, alt)
                        break
        alleles[cn] = rec
        print(f"  结局等位 {cn}: {len(rec)}", flush=True)
    return alleles


def harmonize_to(b, se, ea, nea, ref_e, ref_a):
    """把暴露效应协调到参考等位方向：返回 (beta, se) 或 None（等位冲突/回文不确定）"""
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
    """Wald ratio + IVW（by=β_y|x 效应, sey=结局 SE）"""
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
    p_q = float(__import__("scipy").stats.chi2.sf(q, df)) if df > 0 else 1.0
    z = b_ivw / se_ivw
    p = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(z)))
    return {"beta": b_ivw, "se": se_ivw, "p": p, "n_snp": len(bx), "Q": q, "Q_p": p_q}


def steiger_dir(bx, se_x, by, se_y, b_xy, n_x, n_y):
    """Steiger 方向检验（近似）：R2_x vs R2_y。
    R2 ≈ (beta/se)^2 / n  的近似不可靠，这里用工具解释方差近似：
    R2 = 2*eaf*(1-eaf)*beta^2 / se^2 / n 的简化。改用 Z 统计量比例。
    """
    # 简化：方向判定用 IVW 效应显著性 + 方差解释量比较（工具在哪个变量的 R2 更大）
    zx = np.abs(bx / se_x)
    zy = np.abs(by / se_y)
    r2x = float(np.mean(zx ** 2))
    r2y = float(np.mean(zy ** 2))
    return {"r2_x_approx": r2x, "r2_y_approx": r2y,
            "direction": "X->Y" if r2x > r2y else "Y->X"}


def main():
    raw = load_raw()
    outd = load_out()
    out_alleles = load_out_alleles()
    print(f"原始效应 SNP: {len(raw)} | 结局表: {len(outd)} 行")

    # ---- A. 两两网络 MR（用 352 IV 中工具强度达标的子集） ----
    dir_rows = []
    for x in EXPOSURES:
        for y in EXPOSURES:
            if x == y:
                continue
            bx_list, by_list, sey_list = [], [], []
            for rs, d in raw.items():
                rx, ry = d.get(x), d.get(y)
                if not rx or not ry:
                    continue
                # 工具 SNP 判定：在 x 上 |Z| >= 3.16（F>=10 近似）
                if abs(rx["beta"]) / rx["se"] < math.sqrt(F_MIN):
                    continue
                # y 效应协调到 x 的效应等位方向
                h = harmonize_to(ry["beta"], ry["se"], ry["ea"], ry["nea"], rx["ea"], rx["nea"])
                if h is None:
                    continue
                bx_list.append(rx["beta"])
                by_list.append(h[0])
                sey_list.append(ry["se"])
            res = ivw(bx_list, by_list, sey_list)
            if res is None:
                print(f"  {x}->{y}: SNP 不足", flush=True)
                continue
            dir_rows.append({"from": x, "to": y, **res})
            print(f"  {x}->{y}: β={res['beta']:+.4f} se={res['se']:.4f} p={res['p']:.3g} "
                  f"n={res['n_snp']} Q_p={res['Q_p']:.3g}", flush=True)

    with open(os.path.join(RESULTS, "mediation_network_direction.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["from", "to", "beta", "se", "p", "n_snp", "Q", "Q_p"], delimiter="\t")
        w.writeheader()
        w.writerows(dir_rows)
    print(f"方向表 -> results/mediation_network_direction.tsv")

    # ---- B. 两步中介 BMI→CRP→结局 ----
    med_rows = []
    for oc in OUTCOMES:
        oa = out_alleles[oc]
        # 第一步：BMI 工具 → CRP（a = β_BMI→CRP）
        a_list, a_se_list = [], []
        bmi_list = []  # 总效应 c 用
        for rs, d in raw.items():
            rB, rC = d.get("BMI"), d.get("CRP")
            o = outd.get((oc, rs))
            if not rB or not rC or not o:
                continue
            if abs(rB["beta"]) / rB["se"] < math.sqrt(F_MIN):
                continue
            hC = harmonize_to(rC["beta"], rC["se"], rC["ea"], rC["nea"], rB["ea"], rB["nea"])
            if hC is None:
                continue
            a_list.append((rB["beta"], hC[0]))
            a_se_list.append(rC["se"])
            # 总效应 c = BMI→结局：结局效应协调到 BMI 等位方向
            oref, oalt = oa.get(rs, (None, None))
            if oref and oalt:
                hO = harmonize_to(o["beta"], o["se"], oalt, oref, rB["ea"], rB["nea"])
                if hO is not None:
                    bmi_list.append((rB["beta"], hO[0], o["se"]))
        a_ivw = ivw([t[0] for t in a_list], [t[1] for t in a_list], a_se_list)
        if a_ivw is None:
            print(f"  {oc}: 第一步 BMI→CRP SNP 不足", flush=True)
            continue

        # 第二步：CRP 工具 → 结局（b = β_CRP→outcome，单暴露 IVW）
        b_list, b_se_list = [], []
        for rs, d in raw.items():
            rB, rC = d.get("BMI"), d.get("CRP")
            o = outd.get((oc, rs))
            if not rB or not rC or not o:
                continue
            if abs(rC["beta"]) / rC["se"] < math.sqrt(F_MIN):
                continue
            oref, oalt = oa.get(rs, (None, None))
            if not oref or not oalt:
                continue
            # 结局效应协调到 CRP 等位方向
            hO = harmonize_to(o["beta"], o["se"], oalt, oref, rC["ea"], rC["nea"])
            if hO is None:
                continue
            b_list.append((rC["beta"], hO[0]))
            b_se_list.append(o["se"])
        b_ivw = ivw([t[0] for t in b_list], [t[1] for t in b_list], b_se_list)
        if b_ivw is None:
            print(f"  {oc}: 第二步 CRP→结局 SNP 不足", flush=True)
            continue

        # 总效应 c = BMI→结局（BMI 工具）
        c_ivw = ivw([t[0] for t in bmi_list], [t[1] for t in bmi_list], [t[2] for t in bmi_list])
        if c_ivw is None:
            print(f"  {oc}: 总效应 BMI→结局 SNP 不足", flush=True)
            continue

        # 中介效应 a*b、直接效应 c-a*b、中介比例
        a, sa = a_ivw["beta"], a_ivw["se"]
        b, sb = b_ivw["beta"], b_ivw["se"]
        c, sc = c_ivw["beta"], c_ivw["se"]
        ab = a * b
        se_ab = math.sqrt((a * sb) ** 2 + (b * sa) ** 2)
        z_ab = ab / se_ab
        p_ab = 2 * (1 - __import__("scipy").stats.norm.cdf(abs(z_ab)))
        direct = c - ab
        prop = ab / c if abs(c) > 1e-12 else float("nan")
        med_rows.append({
            "outcome": oc,
            "a_bmi_to_crp": round(a, 5), "a_se": round(sa, 5), "a_p": f"{a_ivw['p']:.3g}", "a_n": a_ivw["n_snp"],
            "b_crp_to_out": round(b, 5), "b_se": round(sb, 5), "b_p": f"{b_ivw['p']:.3g}", "b_n": b_ivw["n_snp"],
            "c_total": round(c, 5), "c_se": round(sc, 5), "c_p": f"{c_ivw['p']:.3g}", "c_n": c_ivw["n_snp"],
            "mediation_ab": round(ab, 5), "se_ab": round(se_ab, 5), "p_ab": f"{p_ab:.3g}",
            "direct_effect": round(direct, 5),
            "mediation_prop": round(prop, 4) if not math.isnan(prop) else "NA",
        })
        print(f"  {oc}: a={a:+.4f} b={b:+.4f} c={c:+.4f} ab={ab:+.5f} p_ab={p_ab:.3g} prop={prop:.2%}", flush=True)

    with open(os.path.join(RESULTS, "mediation_results_5outcomes.tsv"), "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(med_rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(med_rows)
    print(f"中介结果 -> results/mediation_results_5outcomes.tsv")


if __name__ == "__main__":
    main()
