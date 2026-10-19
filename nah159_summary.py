import collections
import pandas as pd

a = pd.read_csv("nah159_audit.csv")
print("cells:", len(a), " columns:", a.col.nunique())
print(a.status.value_counts().to_string())
print()
print("--- SCALE / OTHER STATE cells ---")
for _, r in a[a.status.isin(["SCALE", "OTHER STATE"])].iterrows():
    print(f"  col {r['col']:>2} {r['meth']:<12}{r['fun']:<11} {str(r['mol']):<32} "
          f"excel = {r['excel']:>10.5f}  log = {r['log']:>10.5f}  [{r['status']}] {r['detail']}")
print()
print("--- MISSING by column ---")
for c, g in a[a.status == "MISSING"].groupby("col"):
    print(f"  col {c:>2}: {len(g)} cells")
print()
miss = a[a.status == "MISSING"]
for col in ("S", "T"):
    mols = sorted(miss[miss.col == col].mol)
    print(f"col {col}: {len(mols)} missing ->", ", ".join(mols))
for col in ("AE", "AG"):
    have = sorted(set(a[(a.col == col) & (a.status == "OK")].mol))
    print(f"col {col}: only {len(have)} molecules have a log:", ", ".join(have))
print()
print("--- per-molecule count of non-OK cells (worst 15) ---")
bad = a[a.status != "OK"]
c = collections.Counter(bad.mol)
for m, n in c.most_common(15):
    print(f"  {m:<32} {n}")
