import pandas as pd
import mr_tadf_lib as L
import invest_audit as IA

a = pd.read_csv("invest_audit.csv")
for col, mol in (("I", 45), ("I", 46), ("O", 44), ("Q", 44), ("AE", 37), ("AG", 37)):
    r = a[(a.col == col) & (a.mol == mol)].iloc[0]
    print(f"  {col}{mol:<3} excel = {r['excel']:<10.5f} log = {r['log']:.5f}  {r['status']:<8} "
          f"{r['prov']:<24} {r['detail']}")

have = sorted(IA.SF.get(("XSF-TDA", "ALDA0", "B3LYP"), {}))
print("\nXSF-TDA B3LYP ALDA0 covers ids:", have)
print("ids present in the Peng sheet:", [i for i in range(1, 47) if i != 10])
