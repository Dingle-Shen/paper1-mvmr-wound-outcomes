#!/usr/env python3
"""MVMR 预实验：三暴露 (T2D/BMI/CRP) 联合 IV × 5 创面/感染结局
IV 构建：三暴露各自 p<5e-8 lead SNP（GWAS Catalog 缓存）→ 1Mb 窗口内取三暴露联合最显著者 → 统一 IV 集
数据：ogw.associations 拉 3 研究×并集 SNP 效应；FinnGen 定点取 5 结局
方法：MVMR-IVW（多元 WLS）+ 全局 Q + 条件 F（弱工具）+ 单暴露 Wald（混杂方向参考）
输出: results/mvmr_summary.csv + results/mvmr_iv_effects.csv
"""
import os, sys, json, csv, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import ogw, finngen_fetch

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

EXPOSURES = [
    {"name": "T2D", "ieu": "ebi-a-GCST006867", "gcst": "GCST006867"},
    {"name": "BMI", "ieu": "ebi-a-GCST006368", "gcst": "GCST006368"},
    {"name": "CRP", "ieu": "ebi-a-GCST90029070", "gcst": "GCST90029070"},
]
OUTCOMES = [
    {"ep": "finngen_R9_L12_ULCERLOWLIMB", "name": "下肢溃疡"},
    {"ep": "finngen_R9_L12_INFECT_SKIN", "name": "皮肤感染"},
    {"ep": "finngen_R9_L12_CELLULITIS", "name": "蜂窝织炎"},
    {"ep": "finngen_R9_L12_DECUBITANSULCERANDPRESSURE", "name": "压疮"},
    {"ep": "finngen_R9_L12_HIDRADENITISSUP", "name": "化脓性汗腺炎"},
]


def load_lead_snps(gcst):
    with open(os.path.join(RESULTS, f"gcst_{gcst}_snps.json")) as f:
        return json.load(f)


def fetch_effects(ieu_id, rsids):
    """OpenGWAS associations 分批拉取，返回 {rsid: rec}。"""
    out = {}
    for i in range(0, len(rsids), 60):
        batch = rsids[i:i + 60]
        for attempt in range(3):
            res = ogw.associations(batch, [ieu_id])
            if isinstance(res, dict) and "__error__" in res:
                print(f"    associations 重试{i}: {res.get('__error__')}", flush=True)
                time.sleep(10)
                continue
            break
        else:
            raise RuntimeError(f"associations 失败 {ieu_id} 批{i}")
        for r in res:
            if r.get("id") != ieu_id:
                continue
            for k in ("beta", "se", "p", "eaf", "n"):
                v = r.get(k)
                if isinstance(v, str):
                    try:
                        r[k] = float(v)
                    except ValueError:
                        r[k] = None
            out[r["rsid"]] = r
    return out


def build_union_iv(exp_effects):
    """统一 1Mb 窗口裁剪：每窗取三暴露联合最显著(最小 p)的 SNP。
    exp_effects: {exp_name: {rsid: rec}}；rec 含 p/beta/se/ea/nea/eaf/chr/position(hg19)。
    返回 [rsid]。"""
    # rsid -> 三暴露中有效的记录
    pool = {}
    for en, eff in exp_effects.items():
        for rsid, r in eff.items():
            if r.get("p") is None or r["p"] >= 5e-8:
                continue
            se = r.get("se")
            if not se:
                continue
            f = (r["beta"] / se) ** 2
            if f < 10:
                continue
            chrom = str(r["chr"])
            pos = int(r["position"])
            pool.setdefault((chrom, pos // 1_000_000, rsid), {})[en] = f
    # 每窗口窗内比较：联合最小 p（用任一暴露的最小 F 近似）
    bywin = {}
    for (chrom, w, rsid), fs in pool.items():
        best_f = max(fs.values())
        bywin.setdefault((chrom, w), []).append((best_f, rsid))
    ivs = []
    for win, lst in bywin.items():
        lst.sort(key=lambda t: -t[0])
        ivs.append(lst[0][1])
    return ivs


def harmonize_all(iv_rsids, exp_effects, out_recs, chrom38_map):
    """对每暴露分别等位协调到结局链，返回 MVMR 矩阵。
    保留三暴露与结局均有效且非回文的 SNP。
    chrom38_map: {rsid: (chr38, pos38)}
    返回 (rows, drop_log)"""
    rows = []
    for rsid in iv_rsids:
        c38, p38 = chrom38_map[rsid]
        o = out_recs.get((c38, p38))
        if not o:
            continue
        recs = {}
        ok = True
        for en, eff in exp_effects.items():
            r = eff.get(rsid)
            if not r:
                ok = False
                break
            ea, nea = r["ea"], r["nea"]
            o_alt, o_ref = o["alt"], o["ref"]
            if ea == o_alt:
                b, s = r["beta"], r["se"]
            elif ea == o_ref:
                b, s = -r["beta"], r["se"]
            else:
                c = COMP.get(ea)
                if c == o_alt:
                    b, s = r["beta"], r["se"]
                elif c == o_ref:
                    b, s = -r["beta"], r["se"]
                else:
                    ok = False
                    break
            eaf = r.get("eaf")
            palin = {ea, nea} in ({"A", "T"}, {"C", "G"})
            o_palin = {o_alt, o_ref} in ({"A", "T"}, {"C", "G"})
            eaf_mid = (eaf is not None and 0.3 < eaf < 0.7)
            o_af_mid = 0.3 < o.get("af_alt", 0.5) < 0.7
            if palin and o_palin and eaf_mid and o_af_mid:
                ok = False
                break
            recs[en] = (b, s)
        if not ok:
            continue
        rows.append({
            "rsid": rsid, "chr38": c38, "pos38": p38,
            "beta_out": o["beta"], "se_out": o["sebeta"],
            **{f"beta_{en}": recs[en][0] for en in exp_effects},
            **{f"se_{en}": recs[en][1] for en in exp_effects},
        })
    return rows


def mvmr_ivw(B, y, w):
    """多元 WLS: y ~ B @ theta，权重 w=1/se_out^2。返回 (theta, se_theta, p)。"""
    from scipy import stats
    W = np.diag(w)
    XtWX = B.T @ W @ B
    XtWX_inv = np.linalg.pinv(XtWX)
    theta = XtWX_inv @ (B.T @ W @ y)
    resid = y - B @ theta
    # 模型协方差（固定暴露效应假设下）
    cov = XtWX_inv
    se = np.sqrt(np.diag(cov))
    # Q（多元异质性）
    q = float(resid @ W @ resid)
    df = len(y) - B.shape[1]
    from scipy.stats import chi2
    p_q = float(chi2.sf(q, df))
    p = 2 * stats.norm.sf(np.abs(theta / se))
    return theta, se, p, q, df, p_q


def conditional_f(B, Bse):
    """Sanderson-Windmeijer 式条件 F：每暴露对其他暴露回归后残差方差/其 se²。"""
    out = {}
    ncol = B.shape[1]
    for j in range(ncol):
        others = [i for i in range(ncol) if i != j]
        X_o = np.column_stack([np.ones(len(B)), B[:, others]])
        yj = B[:, j]
        coef, *_ = np.linalg.lstsq(X_o, yj, rcond=None)
        resid = yj - X_o @ coef
        partial = np.sum(resid ** 2)
        # 条件 F = (resid_var_explained) / se_j^2 近似：用 1/se_j^2 的残差和
        f = float(np.sum((resid / Bse[:, j]) ** 2) / (len(B) - ncol))
        out[j] = f
    return out


if __name__ == "__main__":
    t0 = time.time()
    print("[1] 加载三暴露 lead SNP 并构建联合 IV ...", flush=True)
    exp_effects = {}
    for e in EXPOSURES:
        leads = load_lead_snps(e["gcst"])
        rsids = [x["rsid"] for x in leads]
        print(f"    {e['name']}: {len(rsids)} lead SNP，拉效应 ...", flush=True)
        exp_effects[e["name"]] = fetch_effects(e["ieu"], rsids)
        print(f"    {e['name']}: IEU 命中 {len(exp_effects[e['name']])}", flush=True)

    iv_rsids = build_union_iv(exp_effects)
    print(f"    联合 IV(统一 1Mb 裁剪): {len(iv_rsids)} SNP", flush=True)

    # 补拉 IV 集上三暴露效应（并集可能含其他暴露未拉过的 SNP）
    for e in EXPOSURES:
        missing = [r for r in iv_rsids if r not in exp_effects[e["name"]]]
        if missing:
            print(f"    补拉 {e['name']}: {len(missing)} SNP", flush=True)
            exp_effects[e["name"]].update(fetch_effects(e["ieu"], missing))

    # liftover hg19 -> hg38
    print("[2] liftover hg19->hg38 ...", flush=True)
    import local_mr
    chrom38_map = {}
    dropped = 0
    for rsid in iv_rsids:
        r = None
        for e in EXPOSURES:
            if rsid in exp_effects[e["name"]]:
                r = exp_effects[e["name"]][rsid]
                break
        if r is None:
            dropped += 1
            continue
        lp = local_mr.lift(str(r["chr"]), int(r["position"]))
        if lp is None:
            dropped += 1
            continue
        chrom38_map[rsid] = lp
    iv_rsids = [r for r in iv_rsids if r in chrom38_map]
    print(f"    liftover 成功 {len(iv_rsids)}, 丢弃 {dropped}", flush=True)

    # 保存 IV 集与效应快照
    with open(os.path.join(RESULTS, "mvmr_iv_effects.csv"), "w", newline="") as f:
        wcsv = csv.writer(f)
        wcsv.writerow(["rsid", "chr38", "pos38"] +
                      [f"{k}_{v}" for e in EXPOSURES for k, v in (("beta", e["name"]), ("se", e["name"]), ("p", e["name"]), ("eaf", e["name"]), ("ea", e["name"]), ("nea", e["name"]))] +
                      ["note"])
        for rsid in iv_rsids:
            c38, p38 = chrom38_map[rsid]
            row = [rsid, c38, p38]
            for e in EXPOSURES:
                r = exp_effects[e["name"]].get(rsid, {})
                row += [r.get("beta"), r.get("se"), r.get("p"), r.get("eaf"), r.get("ea"), r.get("nea")]
            row.append("")
            wcsv.writerow(row)

    # 逐结局跑 MVMR
    summary = []
    for oc in OUTCOMES:
        print(f"\n=== 结局: {oc['name']} ({oc['ep']}) ===", flush=True)
        fe = finngen_fetch.FinnGenEndpoint(oc["ep"])
        out_recs = fe.fetch_many([chrom38_map[r] for r in iv_rsids])
        print(f"    结局命中 {len(out_recs)}/{len(iv_rsids)}", flush=True)
        rows = harmonize_all(iv_rsids, exp_effects, out_recs, chrom38_map)
        print(f"    协调后有效 {len(rows)}", flush=True)
        if len(rows) < 12:
            print("    !! 有效 SNP 过少，跳过", flush=True)
            continue
        B = np.column_stack([[float(r[f"beta_{e['name']}"]) for r in rows] for e in EXPOSURES])
        Bse = np.column_stack([[float(r[f"se_{e['name']}"]) for r in rows] for e in EXPOSURES])
        y = np.array([float(r["beta_out"]) for r in rows])
        sy = np.array([float(r["se_out"]) for r in rows])
        w = 1 / sy ** 2
        theta, se, p, q, df, p_q = mvmr_ivw(B, y, w)
        cf = conditional_f(B, Bse)
        # 单暴露 Wald（混杂方向参考，不做因果解释）
        uni = {}
        for j, e in enumerate(EXPOSURES):
            wj = w * B[:, j] ** 2
            ratio = float(np.sum(wj * (y / B[:, j])) / np.sum(wj))
            ratio_se = float(1 / np.sqrt(np.sum(wj)))
            from scipy import stats as st
            uni[e["name"]] = (ratio, ratio_se, 2 * st.norm.sf(abs(ratio / ratio_se)))
        print(f"    n_IV={len(rows)}  Q={q:.1f}(df={df}, p={p_q:.3g})", flush=True)
        for j, e in enumerate(EXPOSURES):
            print(f"    MVMR {e['name']:>4s}: β={theta[j]:+.4f} SE={se[j]:.4f} p={p[j]:.3g} | 条件F={cf[j]:.0f} | 单暴露Wald β={uni[e['name']][0]:+.4f} p={uni[e['name']][2]:.3g}", flush=True)
        for j, e in enumerate(EXPOSURES):
            summary.append({
                "outcome": oc["name"], "exposure": e["name"], "n_iv": len(rows),
                "mvmr_beta": round(float(theta[j]), 4), "mvmr_se": round(float(se[j]), 4),
                "mvmr_p": f"{float(p[j]):.3g}", "cond_F": round(float(cf[j]), 0),
                "Q": round(q, 1), "Q_p": f"{p_q:.3g}",
                "uni_beta": round(uni[e["name"]][0], 4), "uni_p": f"{uni[e['name']][2]:.3g}",
            })

    outp = os.path.join(RESULTS, "mvmr_summary.csv")
    if not summary:
        print("\n!! 所有结局均因有效 SNP 不足被跳过，无 summary 输出", flush=True)
        sys.exit(0)
    with open(outp, "w", newline="", encoding="utf-8-sig") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        wcsv.writeheader(); wcsv.writerows(summary)
    print(f"\n已保存: {outp}  (耗时 {time.time()-t0:.0f}s)")
