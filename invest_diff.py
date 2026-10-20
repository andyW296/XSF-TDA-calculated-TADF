"""Diff the Peng sheet before/after the workbook update and list what still disagrees."""
import pandas as pd

old = pd.read_csv("prev_invest_audit.csv")
new = pd.read_csv("invest_audit.csv")
m = old.merge(new, on=["col", "mol"], suffixes=("_old", "_new"))
m["excel_old"] = pd.to_numeric(m["excel_old"], errors="coerce")
m["excel_new"] = pd.to_numeric(m["excel_new"], errors="coerce")
ch = m[(m.excel_old.isna() != m.excel_new.isna())
       | ((m.excel_old - m.excel_new).abs() > 1e-9)
       | (m.status_old != m.status_new)]
print(f"cells changed in the sheet / status: {len(ch)}")
for _, r in ch.iterrows():
    print(f"  col {r['col']:>2} {r['meth_new']:<12}{str(r['var_new'] or '-'):<7}{r['fun_new']:<11} "
          f"{str(r['mol']):<28} excel {r['excel_old']} -> {r['excel_new']}   "
          f"[{r['status_old']} -> {r['status_new']}]")
print()
print("current status counts:", new.status.value_counts().to_dict())
print()
print("--- still disagreeing cells ---")
for _, r in new[new.status.isin(["NO SUPPORT", "OTHER STATE", "SCALE", "MISSING"])].iterrows():
    print(f"  col {r['col']:>2} {r['meth']:<12}{str(r['var'] or '-'):<7}{r['fun']:<11} "
          f"{str(r['mol']):<26} excel = {r['excel']:>10.5f} log = "
          f"{'--' if pd.isna(r['log']) else format(r['log'], '.5f'):>10} "
          f"[{r['status']}] {r['prov']} {r['detail']}")
