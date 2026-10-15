# -*- coding: utf-8 -*-
"""第 1 步核验：用本地 R9 全量文件重跑立项书 15 组 MVMR。

背景（核验方 2026-10-12 意见）：立项书 MVMR（10-07）结局走 finngen_fetch.FinnGenEndpoint
实时拉取；本地 R9 全量文件 10-10 才下载完。单变量 Wald 方向在两版间翻转
（皮肤感染 CRP：立项书 -0.636 vs 新版 +0.185），怀疑立项书结局侧对齐/数据版本有问题。
第 1 步 = 隔离"结局数据版本"变量：IV 与暴露全用立项书同款（mvmr_docA_ivs_rebuilt.tsv +
IEU），只把结局换成本地 R9 全量文件，重跑 15 组 MVMR。

输出：results/mvmr_docA_r9local_summary.csv + 对比 mvmr_docA_summary.csv
"""
import os, sys, gzip, csv, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scipy_shim  # noqa: F401  注册最小 scipy 兼容层（本环境无 scipy）
import numpy as np
import ogw
from scipy import stats as st
from mvmr import mvmr_ivw, conditional_f

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
# R9 结局全量文件：优先用工作区副本，缺失时回退到 E 盘归档
_WORKSPACE_R9 = os.path.join(RESULTS, "finngen_full")
_ARCHIVE_R9 = os.environ.get("MR_FINNGEN_R9_DIR", _WORKSPACE_R9)
FULL = _WORKSPACE_R9 if os.path.isdir(_WORKSPACE_R9) and \
    any(f.endswith(".gz") for f in os.listdir(_WORKSPACE_R9)) else _ARCHIVE_R9
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}

EXPOSURES = [
    {"name": "T2D", "ieu": "ebi-a-GCST006867"},
    {"name": "BMI", "ieu": "ebi-a-GCST006368"},
    {"name": "CRP", "ieu": "ebi-a-GCST90029070"},
]
OUTCOMES = [
    ("L12_INFECT_SKIN", "皮肤感染"),
    ("L12_CELLULITIS", "蜂窝织炎"),
    ("L12_ULCERLOWLIMB", "下肢溃疡"),
    ("L12_CHRONICULCEROFSKIN", "慢性皮肤溃疡"),
    ("L12_DECUBITANSULCERANDPRESSURE", "压疮"),
]
FNAMES = {
    "L12_INFECT_SKIN": "finngen_R9_L12_INFECT_SKIN.gz",
    "L12_CELLULITIS": "finngen_R9_L12_CELLULITIS.gz",
    "L12_ULCERLOWLIMB": "finngen_R9_L12_ULCERLOWLIMB.gz",
    "L12_CHRONICULCEROFSKIN": "finngen_R9_L12_CHRONICULCEROFSKIN.gz",
    "L12_DECUBITANSULCERANDPRESSURE": "finngen_R9_L12_DECUBITANSULCERANDPRESSURE.sumstats.txt.gz",
}
MIN_FULL_SIZE = 700 * 1024 * 1024


def load_ivs():
    """mvmr_docA_ivs_rebuilt.tsv: snp chr pos p"""
    ivs = []
    with open(os.path.join(RESULTS, "mvmr_docA_ivs_rebuilt.tsv"), encoding="utf-8") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) >= 2:
                ivs.append(p[0])
    return ivs


def fetch_effects(ieu_id, rsids):
    out = {}
    for i in range(0, len(rsids), 60):
        batch = rsids[i:i + 60]
        for attempt in range(3):
            res = ogw.associations(batch, [ieu_id])
            if isinstance(res, dict) and "__error__" in res:
                print(f"    retry{i}: {res.get('__error__')}", flush=True)
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


def scan_outcome(path, want):
    """流式扫 R9 文件，按 rsid 匹配，返回 {rsid: (alt, ref, beta, se, pval, af_alt, chrom, pos)}"""
    out = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("#chrom"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 12:
                continue
            rsids, chrom, pos, ref, alt = p[4], p[0], p[1], p[2], p[3]
            pval, beta, se, af_alt = p[6], p[8], p[9], p[10]
            if not rsids:
                continue
            for r in rsids.split(";"):
                r = r.strip()
                if r in want:
                    out[r] = (alt, ref, float(beta), float(se), float(pval),
                              float(af_alt) if af_alt else None, chrom, pos)
                    break
    return out


def harmonize_r9(iv_rsids, exp_effects, out_recs):
    rows = []
    drop = {"no_outcome": 0, "no_exposure": 0, "allele_mismatch": 0, "palindrome": 0}
    for rsid in iv_rsids:
        o = out_recs.get(rsid)
        if not o:
            drop["no_outcome"] += 1
            continue
        alt, ref, b_out, se_out, p_out, af_out, chrom, pos = o
        recs = {}
        ok = True
        for en, eff in exp_effects.items():
            r = eff.get(rsid)
            if not r or r.get("ea") is None:
                drop["no_exposure"] += 1
                ok = False
                break
            ea, nea = r["ea"], r["nea"]
            b, s = r["beta"], r["se"]
            if ea == alt:
                pass
            elif ea == ref:
                b, s = -b, s
            else:
                c = COMP.get(ea)
                if c == alt:
                    pass
                elif c == ref:
                    b, s = -b, s
                else:
                    drop["allele_mismatch"] += 1
                    ok = False
                    break
            eaf = r.get("eaf")
            palin = {ea, nea} in ({"A", "T"}, {"C", "G"})
            o_palin = {alt, ref} in ({"A", "T"}, {"C", "G"})
            eaf_mid = (eaf is not None and 0.3 < eaf < 0.7)
            o_af_mid = (af_out is not None and 0.3 < af_out < 0.7)
            if palin and o_palin and eaf_mid and o_af_mid:
                drop["palindrome"] += 1
                ok = False
                break
            recs[en] = (b, s)
        if not ok:
            continue
        rows.append({
            "rsid": rsid, "chr38": chrom, "pos38": pos,
            "beta_out": b_out, "se_out": se_out, "p_out": p_out,
            **{f"beta_{en}": recs[en][0] for en in exp_effects},
            **{f"se_{en}": recs[en][1] for en in exp_effects},
        })
    return rows, drop


def main():
    t0 = time.time()
    iv_rsids = load_ivs()
    print(f"[1] IV 清单: {len(iv_rsids)} SNP", flush=True)

    exp_effects = {}
    for e in EXPOSURES:
        print(f"[2] 拉 {e['name']} ({e['ieu']}) ...", flush=True)
        exp_effects[e["name"]] = fetch_effects(e["ieu"], iv_rsids)
        print(f"    {e['name']}: IEU 命中 {len(exp_effects[e['name']])}/{len(iv_rsids)}", flush=True)

    summary = []
    for epid, cn in OUTCOMES:
        fname = FNAMES[epid]
        path = os.path.join(FULL, fname)
        if not os.path.exists(path) or os.path.getsize(path) < MIN_FULL_SIZE:
            print(f"[SKIP] {epid} 文件未完成", flush=True)
            continue
        print(f"\n=== {cn} ({epid}) 扫描本地 R9 ...", flush=True)
        ts = time.time()
        out_recs = scan_outcome(path, set(iv_rsids))
        print(f"    结局命中 {len(out_recs)}/{len(iv_rsids)} ({time.time()-ts:.0f}s)", flush=True)

        rows, drop = harmonize_r9(iv_rsids, exp_effects, out_recs)
        print(f"    协调后有效 {len(rows)} | drop={drop}", flush=True)
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
        uni = {}
        for j, e in enumerate(EXPOSURES):
            # 标准单暴露 IVW：以比值自身方差的倒数加权（w·B²）。
            # 原实现点估计仅用结局侧 1/se_out² 加权，与 ratio_se 的加权口径不一致，
            # 会系统性放大 |B| 小的借用位点，产生虚假的方向背离。
            wj = w * B[:, j] ** 2
            ratio = float(np.sum(wj * (y / B[:, j])) / np.sum(wj))
            ratio_se = float(1 / np.sqrt(np.sum(wj)))
            uni[e["name"]] = (ratio, ratio_se, 2 * st.norm.sf(abs(ratio / ratio_se)))
        print(f"    n_IV={len(rows)} Q={q:.1f}(df={df}, p={p_q:.3g})", flush=True)
        for j, e in enumerate(EXPOSURES):
            print(f"    MVMR {e['name']:>4s}: β={theta[j]:+.4f} SE={se[j]:.4f} p={p[j]:.3g} | 条件F={cf[j]:.0f} | 单暴露Wald β={uni[e['name']][0]:+.4f} p={uni[e['name']][2]:.3g}", flush=True)
            summary.append({
                "outcome": cn, "exposure": e["name"], "n_iv": len(rows),
                "mvmr_beta": round(float(theta[j]), 4), "mvmr_se": round(float(se[j]), 4),
                "mvmr_p": f"{float(p[j]):.3g}", "cond_F": round(float(cf[j]), 0),
                "Q": round(q, 1), "Q_p": f"{p_q:.3g}",
                "uni_beta": round(uni[e["name"]][0], 4), "uni_p": f"{uni[e['name']][2]:.3g}",
            })

    outp = os.path.join(RESULTS, "mvmr_docA_r9local_summary.csv")
    if not summary:
        print("!! 无 summary 输出", flush=True)
        sys.exit(0)
    with open(outp, "w", newline="", encoding="utf-8-sig") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(summary[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(summary)
    print(f"\n已保存: {outp} (耗时 {time.time()-t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
