"""Final audit of TADF_all.xlsx sheet 'MR-TADF' against every pyscf log in this repo.

Column recipes (all established from the logs themselves, not assumed):
  C K S AA          TDA      S1 - T1        (tddft.TDA, singlet and triplet)
  D L T AB          TD-DFT   S1 - T1        (tddft.TDDFT, singlet and triplet)
  F H N P V X AD AF USF-TDA  energy of the lowest non-triplet state of the SA_SF_TDA
                             spectrum, measured from the triplet reference state
                             (chosen by <S^2> of each state, not by its position)
  G I O Q W Y AE AG XSF-TDA  the 2nd eigenvalue of the XSF_TDA 'excited energy' array

A cell is classified as
  OK          some log in the repo reproduces the cell with that recipe
  OTHER STATE another eigenvalue of the same run equals the cell
  SCALE       the cell equals an eigenvalue / 1000
  NO SUPPORT  no log in the repo produces this number
  MISSING     the cell is filled but the repo has no data for it
"""
import collections
import pandas as pd

import mr_tadf_lib as L

TD, SF, SFSTATE, inventory = L.collect()
df, names = L.load_excel()
nmol = len(names)
TOL = 5e-4


def runs(cidx, m):
    """[(value, provenance, index, pure)] produced by the shipped logs for that cell."""
    kind, meth, variant, fun = L.COLUMNS[cidx]
    out = []
    if kind == "TD":
        s1 = TD.get((meth, fun, "Singlet"), {})
        t1 = TD.get((meth, fun, "Triplet"), {})
        if s1.get(m) and t1.get(m):
            for e1, n1 in s1[m]:
                for e2, n2 in t1[m]:
                    out.append((e1 - e2, n1 if n1 == n2 else f"{n1} & {n2}", 0, True))
        return out
    d = SF.get((meth, variant, fun), {})
    det = SFSTATE.get((meth, variant, fun), {})
    for k, (nums, src) in enumerate(d.get(m, [])):
        detail = det[m][k][0] if k < len(det.get(m, [])) else []
        if meth == "USF-TDA":
            v, i, pure = L.usf_first_singlet(nums, detail)
            if v is None:                       # no <S^2> printed -> fall back to 3rd state
                v, i, pure = nums[2], 2, False
            out.append((v, src, i, pure))
        else:
            out.append((nums[1], src, 1, False))
    return out


def elements(cidx, m):
    kind, meth, variant, fun = L.COLUMNS[cidx]
    if kind == "TD":
        return []
    d = SF.get((meth, variant, fun), {})
    return [(x, src, i) for nums, src in d.get(m, []) for i, x in enumerate(nums)]


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
        nsrc = len(rr)
        if pd.isna(ev):
            rows.append(dict(base, excel=ev, log=rr[0][0], diff=None, status="EXCEL EMPTY",
                             prov=rr[0][1], detail="", nsrc=nsrc))
            continue
        ev = float(ev)
        ok = [(v, p) for v, p, i, _ in rr if abs(v - ev) <= TOL]
        if ok:
            rows.append(dict(base, excel=ev, log=ok[0][0], diff=ok[0][0] - ev, status="OK",
                             prov=ok[0][1], detail="", nsrc=nsrc))
            continue
        # value equals a different eigenvalue of the same calculation?
        el = [(x, p, i) for x, p, i in elements(cidx, m) if abs(x - ev) <= TOL]
        if el:
            rows.append(dict(base, excel=ev, log=el[0][0], diff=el[0][0] - ev,
                             status="OTHER STATE", prov=el[0][1],
                             detail=f"array[{el[0][2]}]", nsrc=nsrc))
            continue
        sc = [(x, p, i) for x, p, i in elements(cidx, m)
              if x != 0 and abs(x / 1000 - ev) <= 1e-3]
        if sc:
            rows.append(dict(base, excel=ev, log=sc[0][0], diff=sc[0][0] - ev,
                             status="SCALE", prov=sc[0][1],
                             detail=f"array[{sc[0][2]}]/1000", nsrc=nsrc))
            continue
        rows.append(dict(base, excel=ev, log=rr[0][0], diff=rr[0][0] - ev,
                         status="NO SUPPORT", prov=rr[0][1], detail="", nsrc=nsrc))

res = pd.DataFrame(rows)
res.to_csv("mr_tadf_audit.csv", index=False, encoding="utf-8-sig")

R = []
R.append("AUDIT: TADF_all.xlsx sheet 'MR-TADF'  vs  every pyscf log shipped in this repository")
R.append("")
R.append("=" * 120)
R.append("1. LOGS USED")
R.append("=" * 120)
for inv in inventory:
    R.append(f"[{inv['folder']:<8}] {inv['log']:<28} gto.M={inv['sections']:>3} "
             f"TD blocks={inv['td']:>3} SF arrays={inv['sf']:>3}  PySCF={','.join(inv['versions'])}")
R.append("")

R.append("=" * 120)
R.append("2. PER-COLUMN VERDICT")
R.append("=" * 120)
R.append(f"{'col':>3} {'method':<8} {'var':<6} {'functional':<11} {'recipe':<24} "
         f"{'OK':>3} {'state':>5} {'scale':>5} {'nosup':>5} {'miss':>4}")
RECIPE = {"TDA": "S1-T1", "TD-DFT": "S1-T1",
          "USF-TDA": "1st singlet (by <S^2>)", "XSF-TDA": "2nd eigenvalue"}
tot = collections.Counter()
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    c = collections.Counter(res[res["col"] == letter]["status"])
    tot.update(c)
    R.append(f"{letter:>3} {meth:<8} {str(variant or '-'):<6} {fun:<11} {RECIPE[meth]:<24} "
             f"{c['OK']:>3} {c['OTHER STATE']:>5} {c['SCALE']:>5} "
             f"{c['NO SUPPORT']:>5} {c['MISSING']:>4}")
R.append("")
R.append(f"   TOTAL: OK {tot['OK']}, other state {tot['OTHER STATE']}, scale {tot['SCALE']}, "
         f"no support {tot['NO SUPPORT']}, missing {tot['MISSING']}, excel empty "
         f"{tot['EXCEL EMPTY']}")

R.append("")
R.append("=" * 120)
R.append("3. CELLS THAT NO SHIPPED LOG REPRODUCES  (the real problems)")
R.append("=" * 120)
bad = res[res["status"].isin(["NO SUPPORT", "OTHER STATE", "SCALE"])]
for _, r in bad.iterrows():
    R.append(f"  col {r['col']:>2} {r['meth']:<8}{str(r['var'] or '-'):<7}{r['fun']:<11} "
             f"#{r['mol']:>2} {str(r['molecule']):<22} excel = {r['excel']:<10.5f} "
             f"log = {r['log']:<10.5f} diff = {r['diff']:+.5f}")
    R.append(f"        verdict: {r['status']}   from {r['prov']} {r['detail']}")
if len(bad) == 0:
    R.append("  none")

R.append("")
R.append("=" * 120)
R.append("4. CELLS TAKEN FROM A DIFFERENT RUN THAN THE FIRST LOG ON DISK")
R.append("   (the number does exist in the repo, so not an error - but the runs disagree)")
R.append("=" * 120)
cnt = 0
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    for m in range(1, nmol + 1):
        sub = res[(res["col"] == letter) & (res["mol"] == m) & (res["status"] == "OK")]
        if not len(sub) or sub.iloc[0]["nsrc"] < 2:
            continue
        okprov = sub.iloc[0]["prov"]
        rr = runs(cidx, m)
        firstsrc = rr[0][1]
        if okprov != firstsrc:
            others = sorted({f"{p}={v:.5f}" for v, p, i, _ in rr})
            R.append(f"  col {letter:>2} {meth:<8}{str(variant or '-'):<7}{fun:<11} #{m:>2} "
                     f"{str(names[m - 1]):<22} excel uses {okprov}")
            R.append(f"        all runs: {', '.join(others)}")
            cnt += 1
if not cnt:
    R.append("  none")
R.append("")

R.append("=" * 120)
R.append("5. MISSING DATA")
R.append("=" * 120)
R.append("5a. the cell is filled but the repo has no such calculation")
miss = res[res["status"] == "MISSING"]
by = collections.OrderedDict()
for _, r in miss.iterrows():
    by.setdefault((r["col"], r["meth"], r["var"], r["fun"]), []).append(int(r["mol"]))
for (col, meth, var, fun), mols in by.items():
    R.append(f"  col {col:>2} {meth:<8}{str(var or '-'):<7}{fun:<11} {len(mols):>2} cells: {mols}")
if not len(miss):
    R.append("  none")
R.append("")
R.append("5b. entire excel blocks without any log in this repository")
groups = collections.defaultdict(list)
for cidx, g in L.NO_SOURCE.items():
    groups[g].append(L.letter(cidx))
for g, cols in groups.items():
    filled = sum(1 for c in L.NO_SOURCE if L.NO_SOURCE[c] == g
                 for m in range(nmol) if not pd.isna(df.iat[3 + m, c]))
    R.append(f"  {g:<17} {','.join(cols):<28} {filled:>3} filled cells - no log here")
R.append("")
R.append("5c. cells that are empty in the excel although a log exists")
emp = res[res["status"] == "EXCEL EMPTY"]
for _, r in emp.iterrows():
    R.append(f"  col {r['col']} #{r['mol']} {names[int(r['mol'])-1]} log = {r['log']:.5f}")
if not len(emp):
    R.append("  none")

R.append("")
R.append("=" * 120)
R.append("6. DEFINITION NOTE  (why the columns are not all the same quantity)")
R.append("=" * 120)
R.append("  C,K,S,AA (TDA) and D,L,T,AB (TD-DFT) are true gaps  S1 - T1.")
R.append("  F,H,N,P,V,X,AD,AF (USF-TDA) and G,I,O,Q,W,Y,AE,AG (XSF-TDA) are the first singlet of a")
R.append("  spin-flip spectrum measured from the triplet reference state; the spin-flip spectrum")
R.append("  still contains the m_s = 0 component of that triplet (it is not spin complete), so the")
R.append("  number is smaller than a genuine S1 - T1 gap by the residual triplet energy:")
for (kind, meth, variant, fun) in [L.COLUMNS[c] for c in (5, 7, 13, 15, 21, 23, 29, 31)]:
    d = SF.get((meth, variant, fun), {})
    det = SFSTATE.get((meth, variant, fun), {})
    gaps = []
    for m in range(1, nmol + 1):
        for k, (nums, src) in enumerate(d.get(m, [])):
            detail = det[m][k][0] if k < len(det.get(m, [])) else []
            s1 = L.usf_first_singlet(nums, detail)
            tri = [e for n, e, ds in detail
                   if n > 1 and ds is not None and ds > -0.3]     # D<S^2> ~ 0 -> triplet
            if s1[0] is not None and tri:
                gaps.append(tri[0] - s1[0])
                break
    if gaps:
        R.append(f"    USF-TDA {variant:<6}{fun:<11} residual m_s=0 triplet - S1 = "
                 f"{min(gaps):+.3f} .. {max(gaps):+.3f} eV, mean {sum(gaps)/len(gaps):+.3f}"
                 f"   ({len(gaps)} molecules)")
    else:
        R.append(f"    USF-TDA {variant:<6}{fun:<11} (no molecule has both a clean triplet and a "
                 f"singlet in the same run)")

R.append("")
R.append("=" * 120)
R.append("7. THE SAME CALCULATION WAS RUN MORE THAN ONCE WITH DIFFERENT RESULTS")
R.append("=" * 120)
R.append("   (each line = one molecule/variant that appears in several logs; the spread is the")
R.append("    difference between the smallest and largest value the logs give for that cell)")
disagree = []
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    letter = L.letter(cidx)
    for m in range(1, nmol + 1):
        rr = runs(cidx, m)
        if len(rr) < 2:
            continue
        vals = [v for v, p, i, _ in rr]
        spread = max(vals) - min(vals)
        if spread > 1e-3:
            disagree.append((spread, letter, meth, variant, fun, m, names[m - 1],
                             sorted({(p, round(v, 5)) for v, p, i, _ in rr})))
disagree.sort(reverse=True)
for spread, letter, meth, variant, fun, m, nm, vals in disagree:
    R.append(f"  col {letter:>2} #{m:>2} {str(nm):<22} {meth:<8}{str(variant or '-'):<7}{fun:<11}"
             f" spread = {spread:.4f} eV")
    for p, v in vals:
        R.append(f"          {p:<30} {v:.5f}")
if not disagree:
    R.append("  none")


out = "\n".join(R)
open("verify_all_mr_tadf_report.txt", "w", encoding="utf-8").write(out + "\n")
print(out)
