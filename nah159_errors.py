"""Detailed list of every NAH159 cell that disagrees with the logs."""
import pandas as pd
import mr_tadf_lib as L
import nah159_audit as NA

SF, SFSTATE, USTD, df = NA.SF, NA.SFSTATE, NA.USTD, NA.df
row_of = {}
for r in range(3, df.shape[0]):
    v = df.iat[r, 0]
    if pd.notna(v):
        row_of[str(v).strip()] = r

res = pd.read_csv("nah159_audit.csv")
bad = res[res.status.isin(["SCALE", "NO SUPPORT", "OTHER STATE"])]
print(f"共 {len(bad)} 个不一致单元格\n")
for _, r in bad.iterrows():
    letter, mol = r["col"], str(r["mol"])
    cidx = [c for c, _ in L.COLUMNS.items() if L.letter(c) == letter][0]
    kind, meth, variant, fun = L.COLUMNS[cidx]
    ex = float(r["excel"])
    print("=" * 104)
    print(f"列 {letter} = {meth} / {variant} / {fun}   分子 = {mol}")
    print(f"   表格值 = {ex:.6f}   [{r['status']}]   {r['detail']}")
    if kind == "USTD":
        for states, src in USTD.get(fun, {}).get(mol, []):
            print(f"   日志 {src}:")
            for n, e, s in states:
                mark = "  <== 表格值" if abs(e - ex) <= 5e-4 else ""
                print(f"      state {n}: {e:>9.5f} eV   <S^2> = {s}{mark}")
    else:
        for nums, src in SF.get((meth, variant, fun), {}).get(mol, []):
            print(f"   日志 {src}:")
            print("      态: " + " ".join(f"{i+1}:{x:>8.4f}" for i, x in enumerate(nums)))
            hit = [f"态{i+1}({x:.5f})" for i, x in enumerate(nums) if abs(x - ex) <= 5e-4]
            hitN = [f"态{i+1}({x:.5f})" for i, x in enumerate(nums)
                    if x and abs(x / 1000 - ex) <= 2e-5]
            print(f"      表格值在数组里: {hit if hit else '无'}")
            print(f"      表格值×1000 对应: {hitN if hitN else '无'}")
            print(f"      表格值×1000 = {ex*1000:.4f}")
    print()
