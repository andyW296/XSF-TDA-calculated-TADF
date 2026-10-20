import collections
import pandas as pd
import mr_tadf_lib as L

a = pd.read_csv("invest_audit.csv")
CH = [("V", 6), ("N", 22), ("N", 31), ("V", 30), ("X", 30), ("V", 34), ("X", 34),
      ("P", 36), ("V", 37), ("X", 37), ("Y", 37), ("V", 39), ("X", 39),
      ("V", 44), ("X", 44)]
print("the 15 cells the update rewrote, and where the number now comes from:")
for col, mol in CH:
    r = a[(a.col == col) & (a.mol == mol)].iloc[0]
    print(f"  {col}{mol:<3} excel = {r['excel']:<12.8g} log = {r['log']:.5f}  {r['status']:<7}"
          f" {r['prov']:<26} {r['detail']}")

import mr_tadf_lib as L
import invest_audit as IA
for key in [("XSF-TDA", "ALDA0", "B3LYP"), ("XSF-TDA", "MC", "B3LYP")]:
    have = sorted(IA.SF.get(key, {}))
    print(f"{key}: {len(have)} molecules, 45/46 present? "
          f"{45 in have}/{46 in have}; missing ids: {[i for i in range(1, 47) if i != 10 and i not in have]}")
