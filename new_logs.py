import re

PATHS = ("Designed INVEST candidates/USF-TDA/huckel_mol37.log",
         "Designed INVEST candidates/USF-TDA/mol44_pbe0.log3",
         "Designed INVEST candidates/USF-TDA/huckel_mol6.log")
for p in PATHS:
    t = open(p, encoding="utf-8", errors="ignore").read()
    print("=" * 90)
    print(p)
    print("   fun lines :", sorted(set(re.findall(r"funs? = (.+)", t)))[:3])
    print("   XC        :", sorted(set(re.findall(r"XC functionals = (\S+)", t))))
    print("   arrays    :")
    for m in re.finditer(r"^excited energy\s+\[([^\]]*)\]", t, re.M):
        print("      ", " ".join(m.group(1).split())[:104])
    print("   states    :")
    for m in re.finditer(r"^(Excited state|Excited State)\s+(\d+)\s+(-?\d+\.\d+\S*)\s+eV(.*)$",
                         t, re.M):
        print(f"      {m.group(1)} {m.group(2)}: {m.group(3)} eV {m.group(4).strip()[:52]}")
        if int(m.group(2)) >= 3:
            break
