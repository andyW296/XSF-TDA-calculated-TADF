"""Compare TADF_all.xlsx sheet 'NAH159' with the NAH159 log folders.

Same recipes as for MR-TADF / Peng:
  TDA / TD-DFT      S1 - T1
  USF-TD-DFT        states 2,3 -> smaller <S^2>      (only the BHHLYP run exists)
  USF-TDA           states 2,3 -> smaller D<S^2>
  XSF-TDA           2nd excitation energy
  excluded          E, M, AC (USF-TD-DFT / B3LYP, PBE0, CAM-B3LYP)

Molecule identity = the folder name 'NAH159/xyz/NNN_<name>/geom.xyz', matched to the
geometry printed in each gto.M section, and compared with column A of the sheet.
"""
import os
import re
import glob
import collections
import pandas as pd

import mr_tadf_lib as L

BASE = "NAH159"
FOLDERS = [os.path.join(BASE, d) for d in ("TDDFT", "USF-TDA", "XSF-TDA")]
SHEET = "NAH159"
TOL = 5e-4


def build_geo():
    geo2name = {}
    for d in sorted(os.listdir(os.path.join(BASE, "xyz"))):
        p = os.path.join(BASE, "xyz", d, "geom.xyz")
        if not os.path.isfile(p):
            continue
        m = re.match(r"^\d+_(.+)$", d)
        if m:
            geo2name[L.read_xyz(p)] = m.group(1)
    return geo2name


geo = build_geo()
print("geometries:", len(geo))
TD, SF, SFSTATE, USTD, inv = L.collect(folders=FOLDERS, geo2mol=geo)

df = pd.ExcelFile(L.XLSX).parse(SHEET, header=None)
rows = []
for r in range(3, df.shape[0]):
    v = df.iat[r, 0]
    if pd.notna(v):
        rows.append((r, str(v).strip()))
print("molecules in the sheet:", len(rows))

RECIPE = {"TDA": "S1-T1", "TD-DFT": "S1-T1",
          "USF-TD-DFT": "state 2/3 with smaller <S^2>",
          "USF-TDA": "state 2/3 with smaller D<S^2>",
          "XSF-TDA": "2nd eigenvalue"}


def runs(cidx, m):
    kind, meth, variant, fun = L.COLUMNS[cidx]
    out = []
    if kind == "TD":
        s1 = TD.get((meth, fun, "Singlet"), {})
        t1 = TD.get((meth, fun, "Triplet"), {})
        if s1.get(m) and t1.get(m):
            for e1, n1 in s1[m]:
                for e2, n2 in t1[m]:
                    out.append((e1 - e2, n1 if n1 == n2 else f"{n1} & {n2}", 0, "S1-T1"))
        return out
    if kind == "USTD":
        for states, src in USTD.get(fun, {}).get(m, []):
            v, n, s2 = L.pick_state_2_3(states, keep="smaller")
            if v is not None:
                out.append((v, src, n, f"state {n}, <S^2>={s2:.2f}"))
        return out
    d = SF.get((meth, variant, fun), {})
    det = SFSTATE.get((meth, variant, fun), {})
    for k, (nums, src) in enumerate(d.get(m, [])):
        detail = det[m][k][0] if k < len(det.get(m, [])) else []
        if meth == "USF-TDA":
            v, n, s2 = L.pick_state_2_3(detail)
            if v is not None:
                out.append((v, src, n, f"state {n}, D<S^2>={s2:.2f}"))
        elif len(nums) > 1:
            out.append((nums[1], src, 2, "2nd eigenvalue"))
    return out


def elements(cidx, m):
    kind, meth, variant, fun = L.COLUMNS[cidx]
    if kind == "TD":
        return []
    if kind == "USTD":
        return [(e, src, n) for states, src in USTD.get(fun, {}).get(m, [])
                for n, e, s in states]
    d = SF.get((meth, variant, fun), {})
    return [(x, src, i + 1) for nums, src in d.get(m, []) for i, x in enumerate(nums)]


out_rows = []
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    for r, name in rows:
        ev = df.iat[r, cidx]
        rr = runs(cidx, name)
        base = dict(col=letter, meth=meth, var=variant, fun=fun, mol=name)
        if not rr:
            out_rows.append(dict(base, excel=ev, log=None, diff=None,
                                 status="MISSING" if not pd.isna(ev) else "EMPTY-BOTH",
                                 prov="", detail="" ))
            continue
        if pd.isna(ev):
            out_rows.append(dict(base, excel=ev, log=rr[0][0], diff=None,
                                 status="EXCEL EMPTY", prov=rr[0][1], detail=rr[0][3]))
            continue
        ev = float(ev)
        ok = [(v, p, d) for v, p, i, d in rr if abs(v - ev) <= TOL]
        if ok:
            out_rows.append(dict(base, excel=ev, log=ok[0][0], diff=ok[0][0] - ev,
                                 status="OK", prov=ok[0][1], detail=ok[0][2]))
            continue
        el = [(x, p, i) for x, p, i in elements(cidx, name) if abs(x - ev) <= TOL]
        if el:
            out_rows.append(dict(base, excel=ev, log=el[0][0], diff=el[0][0] - ev,
                                 status="OTHER STATE", prov=el[0][1],
                                 detail=f"state {el[0][2]}"))
            continue
        sc = [(x, p, i) for x, p, i in elements(cidx, name)
              if x != 0 and abs(x / 1000 - ev) <= 2e-5]
        if sc:
            out_rows.append(dict(base, excel=ev, log=sc[0][0], diff=sc[0][0] - ev,
                                 status="SCALE", prov=sc[0][1],
                                 detail=f"state {sc[0][2]}/1000"))
            continue
        out_rows.append(dict(base, excel=ev, log=rr[0][0], diff=rr[0][0] - ev,
                             status="NO SUPPORT", prov=rr[0][1], detail=rr[0][3]))

res = pd.DataFrame(out_rows)
res.to_csv("nah159_audit.csv", index=False, encoding="utf-8-sig")

R = []
R.append("AUDIT: TADF_all.xlsx sheet 'NAH159'  vs  the NAH159 log folders")
R.append("rules: TDA/TD-DFT S1-T1 | USF-TD-DFT state2/3 smaller <S^2> | "
         "USF-TDA state2/3 smaller D<S^2> | XSF-TDA 2nd eigenvalue")
R.append("")
R.append("=" * 118)
R.append("1. LOGS USED")
R.append("=" * 118)
for i in inv:
    R.append(f"{i['log']:<30} sections={i['sections']:>4} identified={i['identified']:>4} "
             f"TD={i['td']:>4} SF={i['sf']:>4} USF-TDDFT={i['us']:>3} "
             f"PySCF={','.join(i['versions'])}")
    if i["per"]:
        R.append("     " + ", ".join(f"{k[0]}/{k[1]}:{v}" for k, v in sorted(i["per"].items())))
R.append("")
R.append("=" * 118)
R.append("2. PER-COLUMN VERDICT")
R.append("=" * 118)
R.append(f"{'col':>3} {'method':<12} {'var':<6} {'functional':<11} {'recipe':<30} "
         f"{'OK':>4} {'state':>5} {'scale':>5} {'nosup':>5} {'miss':>4}")
tot = collections.Counter()
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    c = collections.Counter(res[res["col"] == letter]["status"])
    tot.update(c)
    R.append(f"{letter:>3} {meth:<12} {str(variant or '-'):<6} {fun:<11} {RECIPE[meth]:<30} "
             f"{c['OK']:>4} {c['OTHER STATE']:>5} {c['SCALE']:>5} "
             f"{c['NO SUPPORT']:>5} {c['MISSING']:>4}")
R.append("")
R.append(f"   {len(res)} cells in {len(L.COLUMNS)} columns -> OK {tot['OK']}   "
         f"other state {tot['OTHER STATE']}   scale {tot['SCALE']}   "
         f"no support {tot['NO SUPPORT']}   missing {tot['MISSING']}   "
         f"excel empty {tot['EXCEL EMPTY']}")

R.append("")
R.append("=" * 118)
R.append("3. CELLS THAT DO NOT MATCH THE LOGS")
R.append("=" * 118)
bad = res[res["status"].isin(["NO SUPPORT", "OTHER STATE", "SCALE", "MISSING"])]
for _, r in bad.iterrows():
    R.append(f"  col {r['col']:>2} {r['meth']:<12}{str(r['var'] or '-'):<7}{r['fun']:<11} "
             f"{str(r['mol']):<32} excel = "
             f"{'' if pd.isna(r['excel']) else format(r['excel'], '.5f'):>10} "
             f"log = {'' if r['log'] is None else format(r['log'], '.5f'):>10} "
             f"[{r['status']}] {r['prov']} {r['detail']}")
if not len(bad):
    R.append("  none")

R.append("")
R.append("=" * 118)
R.append("4. SUMMARY BY METHOD")
R.append("=" * 118)
g = res.groupby("meth").status.value_counts().unstack(fill_value=0)
R.append(g.to_string())

out = "\n".join(R)
open("nah159_report.txt", "w", encoding="utf-8").write(out + "\n")
print(out)
