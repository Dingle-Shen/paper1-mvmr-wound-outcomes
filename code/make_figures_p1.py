# -*- coding: utf-8 -*-
"""论文统计图绘制（英文 + Arial，300 dpi PNG + PDF）。

All figure text in English, Arial font, 300 dpi PNG + PDF.

Fig 1  Forest plot of the 15 multivariable MR estimates
Fig 2  Single-exposure versus multivariable estimates (direction contrast)

本包只生成论文一的两张图（fig1、fig2）。原工程脚本中的 fig3-fig5
（候选基因筛选漏斗、基因层火山图、功效与工具强度）属分子层探索线，
不属于论文一，未包含在本包的主流程中（函数保留但不调用）。
"""
import os
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch
from matplotlib.lines import Line2D

BASE = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
RES = os.path.join(BASE, "results")
OUT = os.path.join(BASE, "figures")
os.makedirs(OUT, exist_ok=True)

rcParams["font.family"] = "Arial"
rcParams["font.sans-serif"] = ["Arial"]
rcParams["axes.unicode_minus"] = False
rcParams["font.size"] = 9
rcParams["axes.linewidth"] = 0.8
rcParams["savefig.dpi"] = 300
rcParams["figure.dpi"] = 300
rcParams["axes.titlesize"] = 10
rcParams["axes.labelsize"] = 10

C_T2D = "#2471A3"
C_BMI = "#C0392B"
C_CRP = "#7D3C98"
C_NS = "#95A5A6"
C_EXP = {"T2D": C_T2D, "BMI": C_BMI, "CRP": C_CRP}
C_FLIP = "#E67E22"

OUT_ORDER = ["Skin infection", "Cellulitis", "Lower leg ulcer",
             "Chronic skin ulcer", "Pressure ulcer"]
OUT_KEY = {
    "Skin infection": "L12_INFECT_SKIN",
    "Cellulitis": "L12_CELLULITIS",
    "Lower leg ulcer": "L12_ULCERLOWLIMB",
    "Chronic skin ulcer": "L12_CHRONICULCEROFSKIN",
    "Pressure ulcer": "L12_DECUBITANSULCERANDPRESSURE",
}
EXP_ORDER = ["T2D", "BMI", "CRP"]
OUT_CN = {v: k for k, v in OUT_KEY.items()}
OUT_CN.update({"皮肤感染": "Skin infection", "蜂窝织炎": "Cellulitis", "下肢溃疡": "Lower leg ulcer",
               "慢性皮肤溃疡": "Chronic skin ulcer", "压疮": "Pressure ulcer"})


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"), bbox_inches="tight",
                    facecolor="white", dpi=300)
    plt.close(fig)
    print("  saved:", name)


def tidy(ax, left=False, bottom=True):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)


# ==================== Fig 1: Forest plot ====================
def fig1():
    d = pd.read_csv(os.path.join(RES, "mvmr_docA_r9local_summary_corrected.csv"))
    d["out_en"] = d["outcome"].map(OUT_CN)
    d["lo"] = d["mvmr_beta"] - 1.96 * d["mvmr_se"]
    d["hi"] = d["mvmr_beta"] + 1.96 * d["mvmr_se"]

    rows, y, yticks = [], 0.0, []
    groups = []
    for oc in OUT_ORDER:
        sub = d[d["out_en"] == oc]
        if not len(sub):
            continue
        g0 = y
        for ex in EXP_ORDER:
            r = sub[sub["exposure"] == ex]
            if not len(r):
                continue
            r = r.iloc[0]
            rows.append((oc, ex, r["mvmr_beta"], r["lo"], r["hi"], r["mvmr_p"]))
            yticks.append(y)
            y += 1.0
        groups.append((oc, (g0 + y - 1) / 2))
        y += 1.4
    ymax, ymin = y - 1.5, -0.9

    fig, ax = plt.subplots(figsize=(9.2, 5.8))
    for oc, yc in groups:
        idx = [i for i, r in enumerate(rows) if r[0] == oc]
        y0, y1 = yticks[idx[0]] - 0.52, yticks[idx[-1]] + 0.52
        ax.axhspan(y0, y1, color="#F2F4F7", zorder=0)
        ax.text(-1.42, (y0 + y1) / 2, oc, ha="left", va="center",
                fontsize=10.5, fontweight="bold", color="#1A5276", zorder=6)
    for i, (oc, ex, b, lo, hi, p) in enumerate(rows):
        yy = yticks[i]
        ax.text(-0.82, yy, ex, ha="left", va="center", fontsize=9.5,
                color="#2C3E50", zorder=6)
    ax.axvline(0, color="#2C3E50", lw=1.0, zorder=1)

    for i, (oc, ex, b, lo, hi, p) in enumerate(rows):
        yy = yticks[i]
        col = C_EXP[ex] if p < 0.05 else C_NS
        ax.plot([lo, hi], [yy, yy], color=col, lw=1.6, solid_capstyle="round", zorder=3)
        ax.plot([lo, lo], [yy - 0.11, yy + 0.11], color=col, lw=1.2, zorder=3)
        ax.plot([hi, hi], [yy - 0.11, yy + 0.11], color=col, lw=1.2, zorder=3)
        ax.scatter([b], [yy], s=52, marker="s", facecolor=col,
                   edgecolor="white", linewidth=0.8, zorder=4)
        star = "***" if p < 0.001 else ("**" if p < 0.01 else ("*" if p < 0.05 else ""))
        ax.annotate(f"{b:+.3f}", xy=(1.005, yy), xycoords=("axes fraction", "data"),
                    ha="left", va="center", fontsize=8.8, family="Arial",
                    annotation_clip=False)
        ax.annotate(f"{p:.2g}", xy=(1.15, yy), xycoords=("axes fraction", "data"),
                    ha="left", va="center", fontsize=8.8, family="Arial",
                    annotation_clip=False)
        ax.annotate(star, xy=(1.24, yy), xycoords=("axes fraction", "data"),
                    ha="left", va="center", fontsize=9.5, color=col,
                    fontweight="bold", annotation_clip=False)

    ax.set_yticks([])
    ax.set_ylim(ymin, ymax)
    ax.set_xlim(-1.48, 1.10)
    ax.set_xlabel("Effect size $\\beta$ (log OR per unit exposure)", fontsize=10.5)
    for xf, lab in ((1.005, "$\\beta$"), (1.15, "P")):
        ax.annotate(lab, xy=(xf, ymax + 0.30), xycoords=("axes fraction", "data"),
                    ha="left", va="bottom", fontsize=10, fontweight="bold",
                    annotation_clip=False)
    tidy(ax, left=False)
    ax.tick_params(axis="x", labelsize=9)
    handles = [Line2D([], [], marker="s", ls="", ms=7, mfc=C_EXP[e], mec="white",
                      label=e) for e in EXP_ORDER]
    handles.append(Line2D([], [], marker="s", ls="", ms=7, mfc=C_NS, mec="white",
                          label="P $\\geq$ 0.05"))
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.52, -0.045),
               ncol=4, frameon=False, fontsize=9.5, handletextpad=0.5,
               columnspacing=1.8)
    ax.text(-1.48, ymax + 0.30, "* P<0.05   ** P<0.01   *** P<0.001",
            ha="left", va="bottom", fontsize=9)
    save(fig, "fig1_forest_mvmr")


# ==================== Fig 2: Single-exposure vs MVMR ====================
def fig2():
    d = pd.read_csv(os.path.join(RES, "mvmr_docA_r9local_summary_corrected.csv"))
    d["out_en"] = d["outcome"].map(OUT_CN)
    vals = pd.concat([d["uni_beta"], d["mvmr_beta"]])
    xmin = vals.min() - 0.35
    xmax = vals.max() + 0.45
    fig, ax = plt.subplots(figsize=(8.6, 5.8))
    ax.axvline(0, color="#2C3E50", lw=1.0, zorder=1)
    nflip, y = 0, 0.0
    for oc in OUT_ORDER:
        sub = d[d["out_en"] == oc]
        if not len(sub):
            continue
        g0 = y
        for ex in EXP_ORDER:
            r = sub[sub["exposure"] == ex]
            if not len(r):
                continue
            r = r.iloc[0]
            u, m = r["uni_beta"], r["mvmr_beta"]
            flip = (np.sign(u) != np.sign(m)) and abs(u) > 0.05 and abs(m) > 0.02
            if flip:
                ax.annotate("", xy=(m, y), xytext=(u, y),
                            arrowprops=dict(arrowstyle="-|>", color=C_FLIP, lw=1.4,
                                            mutation_scale=12), zorder=5)
                nflip += 1
            ax.scatter([u], [y], s=58, marker="o", facecolor="white",
                       edgecolor="#7F8C8D", linewidth=1.4, zorder=4)
            ax.scatter([m], [y], s=66, marker="s", facecolor=C_EXP[ex],
                       edgecolor="white", linewidth=0.8, zorder=4)
            ax.text(xmin + 0.12, y, ex, ha="left", va="center", fontsize=9.5,
                    color="#2C3E50", zorder=6)
            y += 1
        yc = (g0 + y - 1) / 2
        ax.annotate(oc, xy=(-0.015, yc), xycoords=("axes fraction", "data"),
                    ha="right", va="center", fontsize=10.5, fontweight="bold",
                    color="#1A5276", annotation_clip=False)
        y += 1.3
    ax.set_yticks([])
    ax.set_ylim(-1.0, y - 1.0)
    ax.set_xlim(xmin, xmax)
    ax.set_xlabel("Effect size $\\beta$", fontsize=10.5)
    tidy(ax, left=False)
    h = [Line2D([], [], marker="o", ls="", ms=8, mfc="white", mec="#7F8C8D", mew=1.5,
                label="Single-exposure IVW"),
         Line2D([], [], marker="s", ls="", ms=8, mfc="#5D6D7E", mec="white",
                label="Multivariable MVMR (colour = exposure)")]
    fig.legend(handles=h, loc="lower center", bbox_to_anchor=(0.5, -0.05),
               ncol=2, frameon=False, fontsize=9.5, handletextpad=0.6,
               columnspacing=2.0)
    save(fig, "fig2_direction_flip")


# ==================== Fig 3: Screening funnel ====================
def fig3():
    stages = [
        ("Candidate pool", 885, None),
        ("eQTLGen coverage", 579, 306),
        ("Genome-wide significance", 549, 30),
        ("Instrument strength", 549, 0),
        ("Palindromic check", 477, 72),
        ("LD clumping", 232, 245),
    ]
    cols = ["#5D6D7E", "#5499C7", "#48C9B0", "#F4D03F", "#EB984E", "#E74C3C"]

    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    x0, wmax, step, bh = 3.05, 4.6, 1.05, 0.52
    base = stages[0][1]
    for i, ((name, val, exc), col) in enumerate(zip(stages, cols)):
        w = max(wmax * (val / base) ** 0.5, 1.15)
        y = -i * step
        ax.add_patch(FancyBboxPatch((x0, y - bh / 2), w, bh,
                                    boxstyle="round,pad=0.015,rounding_size=0.06",
                                    facecolor=col, edgecolor="white", lw=1.2, zorder=2))
        ax.text(x0 + w / 2, y, name, ha="center", va="center",
                fontsize=9.5, color="white", fontweight="bold", zorder=3)
        ax.text(x0 + w + 0.18, y, f"{val}", ha="left", va="center",
                fontsize=14, fontweight="bold", color=col, zorder=3, family="Arial")
        if exc:
            ax.text(x0 - 0.18, y, f"excluded {exc}", ha="right", va="center",
                    fontsize=9.5, color="#C0392B", zorder=3, family="Arial")
        if i < len(stages) - 1:
            ax.annotate("", xy=(x0 + w / 2, y - bh / 2 - 0.02),
                        xytext=(x0 + w / 2, y - bh / 2 - 0.01),
                        arrowprops=dict(arrowstyle="-|>", color="#7F8C8D", lw=1.3,
                                        mutation_scale=11))
    ylast = -len(stages) * step
    ax.add_patch(FancyBboxPatch((x0, ylast - bh / 2), 2.3, bh,
                                boxstyle="round,pad=0.015,rounding_size=0.06",
                                facecolor="#154360", edgecolor="white", lw=1.2, zorder=2))
    ax.text(x0 + 1.15, ylast, "Included in main analysis", ha="center", va="center",
            fontsize=9.8, color="white", fontweight="bold", zorder=3)
    ax.text(x0 + 2.48, ylast, "232 genes / 407 tools", ha="left", va="center",
            fontsize=12.5, fontweight="bold", color="#154360", zorder=3, family="Arial")
    ax.text(0.02, -0.13,
            "Thresholds: P < 5x10$^{-8}$ for genome-wide significance, "
            "F > 10 for instrument strength, r$^2$ < 0.001 within 10 Mb for LD clumping",
            transform=ax.transAxes, ha="left", va="center", fontsize=9, color="#566573")
    ax.set_xlim(0, 9.4)
    ax.set_ylim(ylast - 0.75, 0.95)
    ax.axis("off")
    save(fig, "fig3_screening_funnel")


# ==================== Fig 4: Volcano ====================
def fig4():
    d = pd.read_csv(os.path.join(RES, "mr_eqtl_multi_r9_results.tsv"), sep="\t")
    d = d[~d["is_control"].astype(str).str.lower().isin(["true", "1"])]
    d["nlp"] = -np.log10(d["p"].clip(lower=1e-300))
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    ax.axhline(-np.log10(0.05), color="#AAB7B8", lw=1.0, ls="--", zorder=1)
    ax.text(-0.88, -np.log10(0.05) + 0.07, "P = 0.05", ha="left", va="bottom",
            fontsize=9, color="#7F8C8D")
    for oc in OUT_ORDER:
        key = OUT_KEY[oc]
        s = d[d["outcome"] == key]
        if not len(s):
            continue
        ax.scatter(s["beta"], s["nlp"], s=15, alpha=0.6, linewidths=0,
                   label=oc, zorder=2)
    fdr = d[d["q_fdr"] < 0.05]
    if len(fdr):
        ax.scatter(fdr["beta"], fdr["nlp"], s=150, marker="D",
                   facecolor="none", edgecolor="#C0392B", linewidth=2.0, zorder=5)
        for _, r in fdr.iterrows():
            ax.annotate(r["gene"], xy=(r["beta"], r["nlp"]),
                        xytext=(0, 16), textcoords="offset points",
                        ha="center", fontsize=11, fontweight="bold",
                        color="#C0392B", zorder=6)
    ax.axvline(0, color="#D5DBDB", lw=0.8, zorder=1)
    ax.set_xlabel("Effect size $\\beta$ (per SD predicted expression)", fontsize=10.5)
    ax.set_ylabel("-log$_{10}$(P)", fontsize=10.5)
    ax.set_xlim(-0.9, 0.9)
    ax.set_ylim(-0.15, 5.9)
    tidy(ax)
    ax.legend(loc="upper left", frameon=False, fontsize=9.5, ncol=2,
              handletextpad=0.4, columnspacing=1.2, markerscale=1.7)
    ax.text(0.98, 0.03, "Red diamonds: FDR q < 0.05",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=9.5,
            color="#C0392B")
    save(fig, "fig4_volcano_screening")


# ==================== Fig 5: Power ====================
def fig5():
    cases = [("Skin infection", 21308), ("Cellulitis", 6683),
             ("Lower leg ulcer", 3784), ("Chronic skin ulcer", 1840),
             ("Pressure ulcer", 1479)]
    r2 = 0.0338
    beta_true = 0.1035
    from math import erf, sqrt

    def sf(x):
        return 0.5 * (1 - erf(x / sqrt(2)))

    def power(nc, r2v, alpha):
        se = 1.0 / np.sqrt(nc * r2v)
        z = beta_true / se - 1.959964
        return max(0.0, 1 - sf(z))

    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.2),
                             gridspec_kw={"width_ratios": [1.12, 1.0],
                                          "wspace": 0.34})
    ax = axes[0]
    names = [c[0] for c in cases]
    p05 = [power(c[1], r2, 0.05) * 100 for c in cases]
    pbo = [power(c[1], r2, 0.05 / 232) * 100 for c in cases]
    y = np.arange(len(names))
    ax.barh(y - 0.20, p05, height=0.38, color="#5499C7", edgecolor="white",
            linewidth=0.7, label="$\\alpha$ = 0.05", zorder=2)
    ax.barh(y + 0.20, pbo, height=0.38, color="#C0392B", edgecolor="white",
            linewidth=0.7, label="Bonferroni (232 tests)", zorder=2)
    ax.axvline(80, color="#909497", lw=1.0, ls="--", zorder=3)
    ax.text(80, -0.72, "80%", ha="center", va="bottom", fontsize=9, color="#707B7C")
    for i, (a, b) in enumerate(zip(p05, pbo)):
        ax.text(a + 2, i - 0.20, f"{a:.1f}", va="center", fontsize=9, family="Arial")
        ax.text(b + 2, i + 0.20, f"{b:.1f}", va="center", fontsize=9, family="Arial")
    ax.set_yticks(y); ax.set_yticklabels(names, fontsize=9.5)
    ax.set_xlim(0, 100); ax.set_ylim(-0.9, len(names) - 0.3)
    ax.invert_yaxis()
    ax.set_xlabel("Statistical power (%)", fontsize=10)
    ax.set_title("R$^2$ = 3.38% (median of instruments)", fontsize=10, pad=8)
    tidy(ax, left=False)
    ax.tick_params(axis="y", length=0)
    ax.legend(loc="lower right", frameon=False, fontsize=9)

    ax2 = axes[1]
    r2s = np.linspace(0.004, 0.25, 240)
    lc = ["#1A5276", "#5499C7", "#48C9B0", "#EB984E", "#C0392B"]
    for (nm, nc), c in zip(cases, lc):
        ps = [power(nc, r, 0.05) * 100 for r in r2s]
        ax2.plot(r2s * 100, ps, lw=1.8, color=c, label=nm)
    ax2.axhline(80, color="#909497", lw=1.0, ls="--")
    ax2.axvline(r2 * 100, color=C_FLIP, lw=1.2, ls=":")
    ax2.annotate("median R$^2$\n3.38%", xy=(r2 * 100, 82), xytext=(r2 * 100 + 1.0, 40),
                 fontsize=9, color=C_FLIP, ha="left", va="center")
    ax2.set_xlabel("Variance explained R$^2$ (%)", fontsize=10)
    ax2.set_ylabel("Statistical power (%)", fontsize=10)
    ax2.set_ylim(0, 104)
    ax2.set_xlim(0, 25)
    ax2.set_title("$\\alpha$ = 0.05", fontsize=10, pad=8)
    tidy(ax2, left=False)
    ax2.legend(loc="lower right", frameon=False, fontsize=8.8)
    save(fig, "fig5_power_boundary")


if __name__ == "__main__":
    print("Drawing figures (English, Arial) ...")
    # 仅论文一的图：图 1 森林图、图 2 方向对照。
    # 原工程脚本另含 fig3-fig5（候选基因筛选漏斗、火山图、功效边界），
    # 属分子层探索线，不属于论文一，故不在此包中生成。
    for fn in (fig1, fig2):
        fn()
    print("\nOutput:", OUT)
    for f in sorted(os.listdir(OUT)):
        if f.endswith(".png"):
            print("  %-32s %7.0f KB" % (f, os.path.getsize(os.path.join(OUT, f)) / 1024))
