"""Where exactly do the 'no support' numbers come from?  Trace every cell to its log."""
import os
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

print("search space = every 'excited energy [...]' eigenvalue of every spin-flip block and")
print("every 'Excited State N:' line of every TD block, in these files:")
for folder in IA.folders:
    files = [os.path.basename(p) for p in sorted(glob.glob(os.path.join(folder, "*")))
             if os.path.isfile(p) and not p.endswith(".py")]
    print(f"  {folder}  ({len(files)} files)")
    print("     " + ", ".join(files))
print()

res = pd.read_csv("invest_audit.csv")
bad = res[res.status == "NO SUPPORT"]
for _, r in bad.iterrows():
    letter, m = r["col"], int(r["mol"])
    cidx = [c for c, _ in L.COLUMNS.items() if L.letter(c) == letter][0]
    kind, meth, variant, fun = L.COLUMNS[cidx]
    ex = float(r["excel"])
    print("=" * 110)
    print(f"col {letter} = {meth} / {variant} / {fun}   molecule {m}   excel = {ex}")
    d = SF.get((meth, variant, fun), {}).get(m, [])
    det = SFSTATE.get((meth, variant, fun), {}).get(m, [])
    if not d:
        print("   no data")
        continue
    for k, (nums, src) in enumerate(d):
        print(f"   file {src}   states 1-{len(nums)}:")
        print("      E  = " + " ".join(f"{x:>8.4f}" for x in nums))
        if k < len(det) and det[k][0]:
            print("      D<S^2> = " + " ".join(
                (f"{s:>8.2f}" if s is not None else f"{'--':>8}") for _, _, s in det[k][0]))
    print(f"   the audit used: file {r['prov']}, state {r['detail']} -> log = {r['log']:.5f}")
    hits = [f"{src} state {i+1} ({x:.5f})" for nums, src in d
            for i, x in enumerate(nums) if abs(x - ex) <= 5e-4]
    print("   excel value appears as: " + (", ".join(hits) if hits else "NOTHING in these arrays"))
    # other files that contain the same molecule / variant / functional
    others = sorted({src for (nums, src) in d})
    print(f"   files holding this molecule: {others}")
