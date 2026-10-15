# -*- coding: utf-8 -*-
"""最小 scipy 兼容垫片（仅覆盖本项目用到的接口）。

本环境的 python 未安装 scipy，而 mvmr.py / mvmr_docA_r12local.py 依赖
scipy.stats.norm.sf 与 scipy.stats.chi2.sf。此模块提供纯 Python 实现，
使既有脚本无需改动即可运行。

用法：把本文件所在目录加入 sys.path，并 import scipy_shim 触发注册。
或在调用脚本顶部 `import scipy_shim`。
"""
import math
import sys
import types


def _norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def _chi2_sf(x, k):
    """卡方上尾概率（不完全伽马函数）。"""
    x = float(x)
    k = float(k)
    if x <= 0:
        return 1.0
    a, xx = k / 2.0, x / 2.0
    if xx < a + 1:
        term = s = 1.0 / a
        n = 0
        while abs(term) > abs(s) * 1e-14 and n < 2000:
            n += 1
            term *= xx / (a + n)
            s += term
        return max(0.0, min(1.0, 1.0 - s * math.exp(-xx + a * math.log(xx) - math.lgamma(a))))
    tiny = 1e-300
    b0, c0 = xx + 1 - a, 1 / tiny
    d0 = 1 / b0
    h = d0
    for i in range(1, 2000):
        an = -i * (i - a)
        b0 += 2
        d0 = an * d0 + b0
        if abs(d0) < tiny:
            d0 = tiny
        c0 = b0 + an / c0
        if abs(c0) < tiny:
            c0 = tiny
        d0 = 1 / d0
        dl = d0 * c0
        h *= dl
        if abs(dl - 1) < 1e-15:
            break
    return max(0.0, min(1.0, h * math.exp(-xx + a * math.log(xx) - math.lgamma(a))))


def _norm_ppf(p):
    """标准正态分位数（Acklam 近似）。"""
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    plow, phigh = 0.02425, 1 - 0.02425
    if p < plow:
        q = math.sqrt(-2 * math.log(p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > phigh:
        q = math.sqrt(-2 * math.log(1 - p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5
    r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


class _Norm:
    @staticmethod
    def sf(z):
        arr = hasattr(z, "__len__")
        if arr:
            import numpy as _np
            return _np.array([_norm_sf(float(v)) for v in z])
        return _norm_sf(float(z))

    @staticmethod
    def cdf(z):
        arr = hasattr(z, "__len__")
        if arr:
            import numpy as _np
            return _np.array([1.0 - _norm_sf(float(v)) for v in z])
        return 1.0 - _norm_sf(float(z))

    @staticmethod
    def ppf(p):
        return _norm_ppf(float(p))


class _Chi2:
    @staticmethod
    def sf(x, df):
        arr = hasattr(x, "__len__")
        if arr:
            import numpy as _np
            return _np.array([_chi2_sf(float(v), df) for v in x])
        return _chi2_sf(float(x), df)


def install():
    """把 scipy / scipy.stats 注册到 sys.modules。"""
    if "scipy" in sys.modules:
        return
    scipy_mod = types.ModuleType("scipy")
    stats_mod = types.ModuleType("scipy.stats")
    stats_mod.norm = _Norm
    stats_mod.chi2 = _Chi2
    scipy_mod.stats = stats_mod
    sys.modules["scipy"] = scipy_mod
    sys.modules["scipy.stats"] = stats_mod


install()
