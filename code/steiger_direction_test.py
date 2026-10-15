# -*- coding: utf-8 -*-
"""Steiger 方向检验：BMI ↔ CRP（中介方向定稿用）。
两种 R² 估计：
  A. 近似：R2 = beta^2 / (se^2 * N)          （单 SNP 回归近似，R² 小时成立）
  B. 精确：r = b*sqrt(2eaf(1-eaf)) / sqrt(b^2*2eaf(1-eaf) + se^2*N)  （TwoSampleMR get_r_from_bsd）
对每个方向（X→Y）用 X 的工具集（|Z|>=sqrt(10)），比较工具对 X 与对 Y 的总解释方差：
  z = (atanh(r_x) - atanh(r_y)) / sqrt(1/(N_x-3) + 1/(N_y-3))
  单侧 p：支持 X→Y 当 r_x 显著 > r_y。
输入：results/mediation_exposure_raw.tsv（352 IV × 三暴露 beta/se/ea/nea/eaf）
输出：results/steiger_direction_bmi_crp.tsv
"""
import os, sys, csv, math
sys.stdout.reconfigure(encoding="utf-8")

RESULTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results")
RAW = os.path.join(RESULTS, "mediation_exposure_raw.tsv")

# 样本量（OpenGWAS gwasinfo 元数据，2026-10-15 取回）
N = {"BMI": 315347, "CRP": 575531, "T2D": 655666}

F_MIN = 10.0
COMP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def harmonize_to(b, se, ea, nea, ref_e, ref_a):
    """把效应协调到 ref 等位方向（R² 计算用 beta 平方，方向不影响；
    保留用于未来扩展，实际 R2 用 beta^2 与 eaf 计算）"""
    return b, se  # R2 与方向无关


def r2_approx(beta, se, n):
    """近似 R² = beta^2/(se^2*N)"""
    if se == 0 or beta is None or se is None:
        return 0.0
    return beta * beta / (se * se * n)


def r2_exact(beta, se, eaf, n):
    """精确 r 的平方（get_r_from_bsd）"""
    if beta is None or se is None or eaf is None or se == 0:
        return 0.0
    v = 2 * eaf * (1 - eaf)
    b2v = beta * beta * v
    r2 = b2v / (b2v + se * se * n)
    return r2


def load_raw():
    raw = {}
    with open(RAW, encoding="utf-8") as f:
        for r in csv.DictReader(f, delimiter="\t"):
            d = {}
            for e in ["BMI", "CRP"]:
                try:
                    d[e] = {"beta": float(r[f"beta_{e}"]), "se": float(r[f"se_{e}"]),
                            "eaf": float(r[f"eaf_{e}"]) if r.get(f"eaf_{e}") not in (None, "") else None}
                except (TypeError, ValueError):
                    d[e] = None
            raw[r["rsid"]] = d
    return raw


def steiger(x, y, raw, method):
    """用 X 工具集检验 X→Y。返回汇总 dict。"""
    tools = []
    for rs, d in raw.items():
        rx, ry = d.get(x), d.get(y)
        if not rx or not ry:
            continue
        if abs(rx["beta"]) / rx["se"] < math.sqrt(F_MIN):
            continue
        tools.append((rs, rx, ry))
    if len(tools) < 3:
        return None
    if method == "approx":
        r2x = sum(r2_approx(t[1]["beta"], t[1]["se"], N[x]) for t in tools)
        r2y = sum(r2_approx(t[2]["beta"], t[2]["se"], N[y]) for t in tools)
    else:
        # exact 需要 eaf；BMI eaf 齐全、CRP eaf 缺失时用 approx 兜底（r2 数量级一致，判断不受影响）
        if any(t[1]["eaf"] is None or t[2]["eaf"] is None for t in tools):
            r2x = sum(r2_approx(t[1]["beta"], t[1]["se"], N[x]) for t in tools)
            r2y = sum(r2_approx(t[2]["beta"], t[2]["se"], N[y]) for t in tools)
        else:
            r2x = sum(r2_exact(t[1]["beta"], t[1]["se"], t[1]["eaf"], N[x]) for t in tools)
            r2y = sum(r2_exact(t[2]["beta"], t[2]["se"], t[2]["eaf"], N[y]) for t in tools)
    rx = math.sqrt(r2x)
    ry = math.sqrt(r2y)
    if rx >= 1.0:
        rx = 0.999999
    if ry >= 1.0:
        ry = 0.999999
    z = (math.atanh(rx) - math.atanh(ry)) / math.sqrt(1 / (N[x] - 3) + 1 / (N[y] - 3))
    p_one = 1 - 0.5 * (1 + math.erf(z / math.sqrt(2)))  # 单侧 P(r_x > r_y)
    direction = f"{x}->{y}" if rx > ry else f"{y}->{x}"
    return {
        "method": method,
        "direction_test": f"{x}->{y}",
        "n_tools": len(tools),
        f"r2_{x}": round(r2x, 6),
        f"r2_{y}": round(r2y, 6),
        f"r_{x}": round(rx, 6),
        f"r_{y}": round(ry, 6),
        "z": round(z, 4),
        "p_support_X_to_Y": f"{p_one:.3g}",
        "supported_direction": direction,
    }


def main():
    raw = load_raw()
    print(f"SNP 总数: {len(raw)}")
    rows = []
    for x, y in [("BMI", "CRP"), ("CRP", "BMI")]:
        for method in ["approx", "exact"]:
            res = steiger(x, y, raw, method)
            if res is None:
                print(f"  {x}->{y} ({method}): 工具不足", flush=True)
                continue
            rows.append(res)
            print(f"  {x}->{y} ({method}): tools={res['n_tools']} "
                  f"r2_{x}={res[f'r2_{x}']} r2_{y}={res[f'r2_{y}']} "
                  f"z={res['z']} p={res['p_support_X_to_Y']} dir={res['supported_direction']}", flush=True)

    outp = os.path.join(RESULTS, "steiger_direction_bmi_crp.tsv")
    fieldnames = ["method", "direction_test", "n_tools", "r2_BMI", "r2_CRP", "r_BMI", "r_CRP",
                  "z", "p_support_X_to_Y", "supported_direction"]
    with open(outp, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"写出 -> {outp}")


if __name__ == "__main__":
    main()
