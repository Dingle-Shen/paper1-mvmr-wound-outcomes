#!/usr/bin/env python3
"""OpenGWAS 轻量客户端（仅标准库），用于检索 GWAS ID、获取元数据、服务端 harmonise(含 LD clumping)。
依赖：Python 3 + JWT token（在 config.py 中设置）。
"""
import json
import urllib.request
import urllib.error
import urllib.parse
import sys

API = "https://api.opengwas.io/api"
# token 从 config 读取，避免硬编码
try:
    from config import OPENGWAS_TOKEN as TOKEN
except Exception:
    TOKEN = ""


def _hdr():
    h = {"Content-Type": "application/json"}
    if TOKEN:
        h["Authorization"] = f"Bearer {TOKEN}"
    return h


def _post(path, body):
    url = API + path
    data = json.dumps(body).encode()
    req = urllib.request.Request(url, data=data, headers=_hdr(), method="POST")
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return {"__error__": e.code, "__msg__": e.read().decode()[:300]}


def _get(path):
    url = API + path
    req = urllib.request.Request(url, headers=_hdr(), method="GET")
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return json.loads(r.read().decode())
    except urllib.error.HTTPError as e:
        return {"__error__": e.code, "__msg__": e.read().decode()[:300]}


def user():
    return _get("/user")


def search(q, n=30):
    """在 GWAS 目录中搜索 trait，返回精简列表。带 size 上限避免整库拖取卡死。"""
    size = min(n, 50)
    res = _get("/gwasinfo?q=" + urllib.parse.quote(q) + f"&size={size}")
    if isinstance(res, dict) and "__error__" in res:
        return res
    out = []
    for r in res[:n]:
        out.append({
            "id": r.get("id"),
            "trait": r.get("trait"),
            "sample_size": r.get("sample_size"),
            "population": r.get("population"),
            "year": r.get("year"),
            "author": r.get("author"),
            "consortium": r.get("consortium"),
        })
    return out


def gwasinfo(ids):
    if isinstance(ids, str):
        ids = [ids]
    return _post("/gwasinfo", {"id": ids})


def harmonise(exp_id, out_ids, pval=5e-8, r2=0.001, kb=10000, proxies=1):
    """服务端 harmonise + LD clumping。返回按 outcome 组织的协调后数据。"""
    if isinstance(out_ids, str):
        out_ids = [out_ids]
    body = {
        "id": exp_id,
        "pval": pval,
        "r2": r2,
        "kb": kb,
        "proxies": proxies,
        "mode": "forward",
        "outcomes": out_ids,
    }
    return _post("/ld/harmonise", body)


def associations(variantlist, idlist):
    """取给定变异在给定 GWAS 中的效应（非 LD）。variantlist: rsid/chr:pos; idlist: 研究ID。"""
    if isinstance(variantlist, str):
        variantlist = [variantlist]
    if isinstance(idlist, str):
        idlist = [idlist]
    return _post("/associations", {"variant": variantlist, "id": idlist})


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "user"
    if cmd == "user":
        print(json.dumps(user(), ensure_ascii=False, indent=2))
    elif cmd == "search":
        q = sys.argv[2]
        res = search(q)
        print(json.dumps(res, ensure_ascii=False, indent=2))
    elif cmd == "info":
        ids = sys.argv[2].split(",")
        print(json.dumps(gwasinfo(ids), ensure_ascii=False, indent=2))
