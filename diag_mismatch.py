"""Dump the SF arrays behind each disagreeing excel cell and look for the value used."""
import pandas as pd
import mr_tadf_lib as L

TD, SF, SFSTATE, inventory = L.collect()
df, names = L.load_excel()

CASES = [(5, 24), (5, 30), (5, 34),
         (7, 30), (7, 34),
         (6, 23), (6, 25), (6, 30), (6, 34),
         (8, 13), (8, 23), (8, 30), (8, 34),
         (14, 23), (16, 23),
         (21, 35), (23, 35),
         (30, 9), (30, 23), (30, 34),
         (32, 9), (32, 20), (32, 23), (32, 34)]

for cidx, m in CASES:
    kind, meth, variant, fun = L.COLUMNS[cidx]
    ex = df.iat[3 + m - 1, cidx]
    print("=" * 108)
    print(f"col {L.letter(cidx)} = {meth} {variant} {fun} | molecule {m:02d} {names[m-1]} | excel = {ex!r}")
    d = SF.get((meth, variant, fun), {})
    if not d.get(m):
        print("   no log data")
        continue
    for nums, src in d[m]:
        print(f"   {src:<28} " + " ".join(f"{x:>8.5f}" for x in nums))
    arr = d[m][0][0]
    hits = []
    for i, x in enumerate(arr):
        if abs(x - ex) <= 5e-4:
            hits.append(f"array[{i}]")
        if x != 0 and abs(x / 1000 - ex) <= 2e-6:
            hits.append(f"array[{i}]/1000")
        if abs(-x - ex) <= 5e-4:
            hits.append(f"-array[{i}]")
    for i in range(len(arr)):
        for j in range(len(arr)):
            if i > j and abs(arr[i] - arr[j] - ex) <= 5e-4:
                hits.append(f"a[{i}]-a[{j}]")
    # compare with the other variant / other folder for the same molecule
    for tag in ["USF-TDA", "XSF-TDA"]:
        for var in ["ALDA0", "MC"]:
            o = SF.get((tag, var, fun), {}).get(m)
            if o and (tag, var) != (meth, variant):
                print(f"   (cross) {tag} {var:<6}        " +
                      " ".join(f"{x:>8.5f}" for x in o[0][0]))
    print("   -> excel value found as:", ", ".join(hits) if hits else "nothing")
