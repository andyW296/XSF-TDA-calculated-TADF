"""What changed between the two committed versions of TADF_all.xlsx?"""
import pandas as pd
from openpyxl.utils import get_column_letter

old = pd.ExcelFile("old_prev.xlsx")
new = pd.ExcelFile("old_head.xlsx")
print("sheets:", old.sheet_names)
for sh in new.sheet_names:
    a = old.parse(sh, header=None)
    b = new.parse(sh, header=None)
    if a.shape != b.shape:
        print(f"{sh}: shape {a.shape} -> {b.shape}")
        continue
    diff = []
    for r in range(a.shape[0]):
        for c in range(a.shape[1]):
            x, y = a.iat[r, c], b.iat[r, c]
            if (pd.isna(x) != pd.isna(y)) or (not pd.isna(x) and x != y):
                diff.append((get_column_letter(c + 1), r, x, y))
    print(f"{sh}: {len(diff)} changed cells")
    for col, r, x, y in diff:
        name = b.iat[r, 0] if not pd.isna(b.iat[r, 0]) else ""
        print(f"    {col}{r+1:<4} {str(name):<26} {x!r:>22} -> {y!r}")
