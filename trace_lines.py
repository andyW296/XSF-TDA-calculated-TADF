"""Line numbers of the spin-flip arrays in the INVEST logs (for manual checking)."""
import os
import re

BASE = "Designed INVEST candidates"
TARGET = {
    "usf_tda_bhhlyp.log": ["ALDA0", "MC"],
    "alda0_mc_cam_b3lyp.log": ["ALDA0", "MC"],
    "cam_b3lyp_mol6.log": ["ALDA0", "MC"],
    "alda0_b3lyp_pbe0_cam.log": ["ALDA0", "MC"],
    "alda0_pbe0_bhhlyp.log": ["ALDA0"],
}
MOL = re.compile(r"^(\d+)_optimized\.xyz")
VAR = re.compile(r"(ALDA0|Multicollinear)")


def find(name):
    for sub in os.listdir(BASE):
        p = os.path.join(BASE, sub, name)
        if os.path.isfile(p):
            return p
    return None


for name, variants in TARGET.items():
    path = find(name)
    if not path:
        continue
    print("=" * 100)
    print(path)
    curvar, curmol = None, None
    for i, line in enumerate(open(path, encoding="utf-8", errors="ignore"), 1):
        s = line.strip()
        m = MOL.match(s)
        if m:
            curmol = int(m.group(1))
            continue
        if not s.startswith("method=0") and "Multicollinear" in s and "ALDA0" not in s:
            curvar = "MC"
        elif "ALDA0" in s and "Multicollinear" not in s:
            curvar = "ALDA0"
        if s.startswith("excited energy"):
            if curmol in (6, 30, 34, 37, 39, 44) or (curmol in (22, 31, 36) and curvar):
                print(f"   line {i:>6}  molecule {curmol:>2}  {curvar}  | {s[:78]}")
