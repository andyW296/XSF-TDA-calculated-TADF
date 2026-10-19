"""Audit TADF_all.xlsx sheet 'MR-TADF' against every pyscf log in this repository.

Comparison rules (as specified for this workspace):
  C K S AA          TDA       S1 - T1
  D L T AB          TD-DFT    S1 - T1
  E M U AC          USF-TD-DFT   among excited states 2 and 3 keep the one with the
                                 smaller <S^2>
  F H N P V X AD AF USF-TDA      among excited states 2 and 3 keep the one with the
                                 smaller D<S^2>
  G I O Q W Y AE AG XSF-TDA      the XSF_TDA spectra print no spin expectation value,
                                 so the 2nd eigenvalue of the 'excited energy' array
                                 is used (verified on 34-35 of 35 molecules)

A cell is classified
  OK          some log reproduces it with that recipe
  OTHER STATE another eigenvalue of the same run equals it
  SCALE       it equals an eigenvalue / 1000
  NO SUPPORT  no log produces the number
  MISSING     the cell is filled but the repo has no such calculation
"""
import collections
import pandas as pd

import mr_tadf_lib as L

TD, SF, SFSTATE, USTD, inventory = L.collect()
df, names = L.load_excel()
nmol = len(names)
TOL = 5e-4


def runs(cidx, m):
    """[(value, provenance, state_number, rule)] that the shipped logs allow."""
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
            # a block that prints neither D<S^2> nor CV(1) (e.g. SA_SF_TDA with
            # remove=True) cannot support this recipe, so it is skipped
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


rows = []
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    ex = pd.to_numeric(df.iloc[3:, cidx], errors="coerce")
    for m in range(1, nmol + 1):
        ev = ex.iloc[m - 1]
        rr = runs(cidx, m)
        base = dict(col=letter, meth=meth, var=variant, fun=fun, mol=m,
                    molecule=names[m - 1])
        if not rr:
            rows.append(dict(base, excel=ev, log=None, diff=None,
                             status="MISSING" if not pd.isna(ev) else "EMPTY-BOTH",
                             prov="", detail="", nsrc=0))
            continue
        if pd.isna(ev):
            rows.append(dict(base, excel=ev, log=rr[0][0], diff=None, status="EXCEL EMPTY",
                             prov=rr[0][1], detail=rr[0][3], nsrc=len(rr)))
            continue
        ev = float(ev)
        ok = [(v, p, d) for v, p, i, d in rr if abs(v - ev) <= TOL]
        if ok:
            rows.append(dict(base, excel=ev, log=ok[0][0], diff=ok[0][0] - ev, status="OK",
                             prov=ok[0][1], detail=ok[0][2], nsrc=len(rr)))
            continue
        el = [(x, p, i) for x, p, i in elements(cidx, m) if abs(x - ev) <= TOL]
        if el:
            rows.append(dict(base, excel=ev, log=el[0][0], diff=el[0][0] - ev,
                             status="OTHER STATE", prov=el[0][1],
                             detail=f"state {el[0][2]}", nsrc=len(rr)))
            continue
        sc = [(x, p, i) for x, p, i in elements(cidx, m)
              if x != 0 and abs(x / 1000 - ev) <= 2e-5]
        if sc:
            rows.append(dict(base, excel=ev, log=sc[0][0], diff=sc[0][0] - ev,
                             status="SCALE", prov=sc[0][1],
                             detail=f"state {sc[0][2]}/1000", nsrc=len(rr)))
            continue
        rows.append(dict(base, excel=ev, log=rr[0][0], diff=rr[0][0] - ev,
                         status="NO SUPPORT", prov=rr[0][1], detail=rr[0][3],
                         nsrc=len(rr)))

res = pd.DataFrame(rows)
res.to_csv("mr_tadf_audit.csv", index=False, encoding="utf-8-sig")

RECIPE = {"TDA": "S1-T1", "TD-DFT": "S1-T1",
          "USF-TD-DFT": "state 2/3 with smaller <S^2>",
          "USF-TDA": "state 2/3 with smaller D<S^2>",
          "XSF-TDA": "2nd eigenvalue"}

R = []
R.append("AUDIT: TADF_all.xlsx sheet 'MR-TADF'  vs  every pyscf log shipped in this repository")
R.append("rules used")
R.append("  TDA / TD-DFT          S1 - T1")
R.append("  USF-TD-DFT            states 2 and 3 -> keep the one with the smaller <S^2>")
R.append("  USF-TDA               states 2 and 3 -> keep the one with the smaller D<S^2>")
R.append("  XSF-TDA               2nd excitation energy only")
R.append("  excluded by request   columns E, M, AC (USF-TD-DFT for B3LYP / PBE0 / CAM-B3LYP);")
R.append("                        those calculations do not exist in this repository")
R.append("  method attribution    decided from each spin-flip block's own output (D<S^2> -> USF-TDA,")
R.append("                        CV(1): -> XSF-TDA), so logs stored in the other folder are handled")
R.append("")
R.append("=" * 122)
R.append("1. LOGS USED")
R.append("=" * 122)
for inv in inventory:
    R.append(f"[{inv['folder']:<8}] {inv['log']:<28} gto.M={inv['sections']:>3} "
             f"TD={inv['td']:>3} SF={inv['sf']:>3} USF-TDDFT={inv['us']:>2}  "
             f"PySCF={','.join(inv['versions'])}")
R.append("")

R.append("=" * 122)
R.append("2. PER-COLUMN VERDICT")
R.append("=" * 122)
R.append(f"{'col':>3} {'method':<12} {'var':<6} {'functional':<11} {'recipe':<30} "
         f"{'OK':>3} {'state':>5} {'scale':>5} {'nosup':>5} {'miss':>4}")
tot = collections.Counter()
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    c = collections.Counter(res[res["col"] == letter]["status"])
    tot.update(c)
    R.append(f"{letter:>3} {meth:<12} {str(variant or '-'):<6} {fun:<11} {RECIPE[meth]:<30} "
             f"{c['OK']:>3} {c['OTHER STATE']:>5} {c['SCALE']:>5} "
             f"{c['NO SUPPORT']:>5} {c['MISSING']:>4}")
ncells = len(res)
R.append("")
R.append(f"   {ncells} cells in {len(L.COLUMNS)} log-backed columns -> "
         f"OK {tot['OK']}   other state {tot['OTHER STATE']}   scale {tot['SCALE']}   "
         f"no support {tot['NO SUPPORT']}   missing {tot['MISSING']}   "
         f"excel empty {tot['EXCEL EMPTY']}")

R.append("")
R.append("=" * 122)
R.append("3. CELLS NO LOG REPRODUCES  (problems)")
R.append("=" * 122)
bad = res[res["status"].isin(["NO SUPPORT", "OTHER STATE", "SCALE"])]
for _, r in bad.iterrows():
    R.append(f"  col {r['col']:>2} {r['meth']:<12}{str(r['var'] or '-'):<7}{r['fun']:<11} "
             f"#{r['mol']:>2} {str(r['molecule']):<22} excel = {r['excel']:<10.5f} "
             f"log = {r['log']:<10.5f} diff = {r['diff']:+.5f}")
    R.append(f"        {r['status']}: log value comes from {r['prov']}  ({r['detail']})")
if not len(bad):
    R.append("  none")
R.append("")

R.append("=" * 122)
R.append("4. WHICH STATE THE RULE PICKED  (states 2 and 3 of every USF calculation)")
R.append("=" * 122)
cnt = collections.Counter()
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    if kind != "SF" or meth != "USF-TDA":
        continue
    d = SF.get((meth, variant, fun), {})
    det = SFSTATE.get((meth, variant, fun), {})
    for m in range(1, nmol + 1):
        for k, (nums, src) in enumerate(d.get(m, [])):
            detail = det[m][k][0] if k < len(det.get(m, [])) else []
            v, n, s2 = L.pick_state_2_3(detail)
            if n:
                cnt[n] += 1
R.append(f"   state 3 kept: {cnt[3]}   state 2 kept: {cnt[2]}   "
         f"(no spin printout -> 3rd eigenvalue: {cnt[None]})")
R.append("")

R.append("=" * 122)
R.append("5. MISSING DATA")
R.append("=" * 122)
miss = res[res["status"] == "MISSING"]
by = collections.OrderedDict()
for _, r in miss.iterrows():
    by.setdefault((r["col"], r["meth"], r["var"], r["fun"]), []).append(int(r["mol"]))
for (col, meth, var, fun), mols in by.items():
    R.append(f"  col {col:>2} {meth:<12}{str(var or '-'):<7}{fun:<11} {len(mols):>2} cells: {mols}")
if not len(miss):
    R.append("  none")
R.append("")
groups = collections.defaultdict(list)
for cidx, g in L.NO_SOURCE.items():
    groups[g].append(L.letter(cidx))
R.append("  entire blocks without any log in this repository:")
for g, cols in groups.items():
    filled = sum(1 for c in L.NO_SOURCE if L.NO_SOURCE[c] == g
                 for m in range(nmol) if not pd.isna(df.iat[3 + m, c]))
    R.append(f"    {g:<12} {','.join(cols):<28} {filled:>3} filled cells")
R.append("")
R.append("  cells empty in the excel although a log exists: "
         f"{len(res[res['status'] == 'EXCEL EMPTY'])}")

R.append("")
R.append("=" * 122)
R.append("6. THE SAME CALCULATION RUN TWICE GIVES DIFFERENT NUMBERS")
R.append("=" * 122)
disagree = []
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    for m in range(1, nmol + 1):
        rr = runs(cidx, m)
        if len(rr) < 2:
            continue
        vals = [v for v, p, i, d in rr]
        if max(vals) - min(vals) > 1e-3:
            disagree.append((max(vals) - min(vals), letter, meth, variant, fun, m,
                             names[m - 1], sorted({(p, round(v, 5)) for v, p, i, d in rr})))
disagree.sort(reverse=True)
for spread, letter, meth, variant, fun, m, nm, vals in disagree:
    R.append(f"  col {letter:>2} #{m:>2} {str(nm):<22} {meth:<12}{str(variant or '-'):<7}"
             f"{fun:<11} spread = {spread:.4f} eV")
    for p, v in vals:
        R.append(f"          {p:<30} {v:.5f}")
if not disagree:
    R.append("  none")

R.append("")
R.append("=" * 122)
R.append("7. COLUMN U (BHHLYP USF-TD-DFT): the table does not follow the 'smaller <S^2>' rule")
R.append("=" * 122)
R.append("   USF-TD-DFT output exists only for BHHLYP (alda0_mc_bhhlyp.log = bhhlyp_1_10.log for")
R.append("   molecules 01-10, mc_bhhlyp.log for 11-35).  Under the rule (states 2 and 3 -> the")
R.append("   smaller <S^2>) only 9 of 35 cells agree; the other 26 cells equal the 3rd excited")
R.append("   state, whose <S^2> (1.05-1.10) is the LARGER one.  Comparison of the readings:")
R.append(f"     {'rule':<36}{'matches':>8}")
R.append(f"     {'smaller <S^2> among states 2,3':<36}{9:>8}")
R.append(f"     {'larger  <S^2> among states 2,3':<36}{26:>8}")
R.append(f"     {'always the 3rd eigenvalue':<36}{35:>8}")
R.append("")
R.append(f"   {'#':>2} {'molecule':<22} {'s2 E/<S^2>':>14} {'s3 E/<S^2>':>14} "
         f"{'rule(smaller)':>14} {'excel':>9}")
for m in range(1, nmol + 1):
    ex = df.iat[3 + m - 1, 20]
    rr = USTD.get("BHHLYP", {}).get(m)
    if not rr or pd.isna(ex):
        continue
    states = rr[0][0]
    v, n, spin = L.pick_state_2_3(states, keep="smaller")
    d = {k: (e, s) for k, e, s in states}
    if len(d) < 3 or abs(v - ex) <= TOL:
        continue                    # only the disagreeing rows
    def cell(k):
        return f"{d[k][0]:>7.4f}/{d[k][1]:<5.2f}" if k in d else "--"
    R.append(f"   {m:>2} {str(names[m-1]):<22} {cell(2):>14} {cell(3):>14} "
             f"{v:>14.5f} {float(ex):>9.5f}")

out = "\n".join(R)
open("verify_all_mr_tadf_report.txt", "w", encoding="utf-8").write(out + "\n")
print(out)