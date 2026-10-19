"""Full provenance of every remaining disagreeing Peng cell: file + line + array."""
import os
import re
import glob
import pandas as pd
import mr_tadf_lib as L
import invest_audit as IA

SF, SFSTATE, df = IA.SF, IA.SFSTATE, IA.df
row_of = {}
for r in range(3, df.shape[0]):
    v = df.iat[r, 0]
    if pd.notna(v):
        row_of[int(v)] = r

# ---- line number of every spin-flip array, per file / molecule / variant ---------
MOL = re.compile(r"(\d+)_optimized\.xyz")
lines_of = {}
for folder in IA.folders:
    for path in sorted(glob.glob(os.path.join(folder, "*"))):
        if not os.path.isfile(path) or path.endswith(".py"):
            continue
        name = os.path.basename(path)
        curvar, curmol = None, None
        for i, line in enumerate(open(path, encoding="utf-8", errors="ignore"), 1):
            s = line.strip()
            m = MOL.search(s)
            if m:
                curmol = int(m.group(1))
            if "Multicollinear" in s and "ALDA0" not in s and not s.startswith("method=0"):
                curvar = "MC"
            elif "ALDA0" in s and "Multicollinear" not in s:
                curvar = "ALDA0"
            if s.startswith("excited energy") and curmol is not None:
                lines_of.setdefault((name, curmol, curvar), []).append(i)

res = pd.read_csv("invest_audit.csv")
bad = res[res.status.isin(["NO SUPPORT", "OTHER STATE"])]
out = []
out.append("provenance of the disagreeing Peng cells")
out.append("")
for _, r in bad.iterrows():
    letter, m = r["col"], int(r["mol"])
    cidx = [c for c, _ in L.COLUMNS.items() if L.letter(c) == letter][0]
    kind, meth, variant, fun = L.COLUMNS[cidx]
    ex = float(r["excel"])
    d = SF.get((meth, variant, fun), {}).get(m, [])
    det = SFSTATE.get((meth, variant, fun), {}).get(m, [])
    out.append("=" * 108)
    out.append(f"col {letter} = {meth} / {variant} / {fun}   molecule {m}   "
               f"excel = {ex:.5f}   [{r['status']}]")
    for k, (nums, src) in enumerate(d):
        ln = lines_of.get((src, m, variant))
        where = f"lines {ln}" if ln else "line ?"
        out.append(f"   {src}  ({where})")
        out.append("      E      : " + " ".join(f"{x:>8.4f}" for x in nums))
        if k < len(det) and det[k][0]:
            ss = [s for _, _, s in det[k][0]]
            out.append("      D<S^2> : " + " ".join(
                (f"{s:>8.2f}" if s is not None else f"{'--':>8}") for s in ss))
    hits = [f"{src} state {i+1} = {x:.5f}" for nums, src in d
            for i, x in enumerate(nums) if abs(x - ex) <= 5e-4]
    out.append(f"   excel {ex:.5f} found as: " +
               (", ".join(hits) if hits else "NO state of these arrays"))
    sel = IA.runs(cidx, m)
    if sel:
        v, prov, n, rule = sel[0]
        out.append(f"   rule picked  : {v:.5f}  ({rule}, file {prov})")
out.append("")
out.append("summary of the rule for each cell:")
out.append(f"  {'col':>3} {'mol':>3} {'excel':>9} {'rule value':>10} {'diff':>8}  file:line")
for _, r in bad.iterrows():
    letter, m = r["col"], int(r["mol"])
    cidx = [c for c, _ in L.COLUMNS.items() if L.letter(c) == letter][0]
    kind, meth, variant, fun = L.COLUMNS[cidx]
    sel = IA.runs(cidx, m)
    if not sel:
        continue
    v, prov, n, rule = sel[0]
    ln = lines_of.get((prov, m, variant))
    out.append(f"  {letter:>3} {m:>3} {float(r['excel']):>9.5f} {v:>10.5f} "
               f"{v-float(r['excel']):>+8.5f}  {prov}:{ln}")

text = "\n".join(out)
open("invest_provenance.txt", "w", encoding="utf-8").write(text + "\n")
print(text)
