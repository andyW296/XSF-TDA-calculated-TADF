"""For each disagreeing excel cell, search every number printed in every MR-TADF log."""
import os
import re
import glob

TARGETS = [
    ("G", 25, 0.00191), ("W", 34, 0.24847),
    ("AE", 9, 0.33757), ("AE", 34, 0.30521),
    ("AG", 9, 0.54019), ("AG", 20, 0.36530), ("AG", 34, 0.51502),
    # values the state-2/3 rule gives, for comparison
    ("U", 3, 1.00101),
]
NUM = re.compile(r"-?\d+\.\d+")
files = []
for d in ("MR-TADF/TDDFT", "MR-TADF/USF-TDA", "MR-TADF/XSF-TDA"):
    files += sorted(glob.glob(os.path.join(d, "*")))

cache = {}
for f in files:
    if os.path.isfile(f):
        cache[os.path.basename(f)] = open(f, encoding="utf-8", errors="ignore").read()

for col, mol, val in TARGETS:
    print("=" * 96)
    print(f"excel {col}{mol} = {val}")
    hits = []
    for name, txt in cache.items():
        for line in txt.splitlines():
            nums = [float(x) for x in NUM.findall(line)]
            for x in nums:
                if abs(x - val) <= 1e-4:
                    hits.append((name, x, line.strip()[:110]))
                    break
            if len(hits) > 6:
                break
        if len(hits) > 6:
            break
    if hits:
        for name, x, line in hits[:6]:
            print(f"   {name:<26} {x:<10} | {line}")
    else:
        # also try x1000 and x/1000 in case of a scaling slip
        alt = []
        for name, txt in cache.items():
            for line in txt.splitlines():
                nums = [float(x) for x in NUM.findall(line)]
                for x in nums:
                    if abs(x / 1000 - val) <= 1e-3 or abs(x * 1000 - val) <= 1e-3:
                        alt.append((name, x, line.strip()[:110]))
                        break
        if alt:
            print("   no exact match; scaling candidates:")
            for name, x, line in alt[:5]:
                print(f"   {name:<26} {x:<10} | {line}")
        else:
            print("   NOT FOUND anywhere in the MR-TADF logs")
