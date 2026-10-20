import re
import pandas as pd
import mr_tadf_lib as L
import invest_audit as IA

SF, SFSTATE, df = IA.SF, IA.SFSTATE, IA.df
row_of = {}
for r in range(3, df.shape[0]):
    v = df.iat[r, 0]
    if pd.notna(v):
        row_of[int(v)] = r

print("--- imports of the two new files ---")
for p in ("Designed INVEST candidates/USF-TDA/huckel_mol37.log",
          "Designed INVEST candidates/USF-TDA/mol44_pbe0.log3",
          "Designed INVEST candidates/USF-TDA/huckel_mol6.log"):
    t = open(p, encoding="utf-8", errors="ignore").read()
    print(f"{p}:  SF_XTDA={'SF_XTDA' in t}  XSF_TDA={'XSF_TDA' in t}  "
          f"has D<S^2>={'D<S^2>' in t}  has CV(={bool(re.search(r'^CV[(]', t, re.M))}")

print()
print("--- where did those arrays land? ---")
for meth in ("USF-TDA", "XSF-TDA"):
    for var in ("ALDA0", "MC"):
        for fun in ("CAM-B3LYP", "PBE0", "BHHLYP", "B3LYP"):
            d = SF.get((meth, var, fun), {})
            for mol in (6, 37, 44):
                for nums, src in d.get(mol, []):
                    if "huckel" in src or "log3" in src:
                        print(f"  [{meth} {var} {fun}] mol {mol}: {src} -> "
                              + " ".join(f"{x:.5f}" for x in nums))

print()
print("--- sheet vs audit for the four cells ---")
for col, mol in (("O", 44), ("Q", 44), ("AE", 37), ("AG", 37),
                 ("N", 44), ("P", 44), ("AD", 37), ("AF", 37)):
    cidx = [c for c, _ in L.COLUMNS.items() if L.letter(c) == col][0]
    kind, meth, variant, fun = L.COLUMNS[cidx]
    ex = df.iat[row_of[mol], cidx]
    rr = IA.runs(cidx, mol)
    print(f"  col {col}{mol} ({meth} {variant} {fun}): excel = {ex}  rule -> "
          + (f"{rr[0][0]:.5f} ({rr[0][1]})" if rr else "no data"))
