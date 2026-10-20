import pandas as pd
import mr_tadf_lib as L

df = pd.ExcelFile(L.XLSX).parse("Peng", header=None)
name_of = {}
for r in range(3, df.shape[0]):
    v = df.iat[r, 0]
    if pd.notna(v):
        name_of[int(v)] = r
a = pd.read_csv("invest_audit.csv")
bad = a[a.status.isin(["NO SUPPORT", "MISSING", "OTHER STATE", "SCALE"])]
print(f"{'col':>3} {'mol':>3} {'excel':>10} {'log':>10}  status      detail")
for _, r in bad.iterrows():
    m = int(r["mol"])
    print(f"{r['col']:>3} {m:>3} {r['excel']:>10.5f} "
          f"{'--' if pd.isna(r['log']) else format(r['log'], '.5f'):>10}  "
          f"{r['status']:<11} {r['prov']} {r['detail']}")
print()
print("column A molecule ids present:", sorted(name_of)[:5], "...", sorted(name_of)[-4:])
for m in sorted({int(x) for x in bad.mol}):
    print(f"  id {m} -> row {name_of[m]+1} (col A value {df.iat[name_of[m],0]})")
