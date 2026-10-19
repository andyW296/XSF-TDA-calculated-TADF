"""Survey MR-TADF/USF-TDA and MR-TADF/XSF-TDA logs."""
import os, re, glob

for D in ["MR-TADF/USF-TDA", "MR-TADF/XSF-TDA"]:
    print("=" * 100)
    print(D)
    print("=" * 100)
    for path in sorted(glob.glob(os.path.join(D, "*"))):
        if not os.path.isfile(path):
            continue
        txt = open(path, encoding="utf-8", errors="ignore").read()
        lines = txt.splitlines()
        end = next((i for i, l in enumerate(lines) if "input file end" in l), 200)
        head = lines[:end]
        funs = [l.strip() for l in head if re.match(r"^\s*funs\s*=", l)]
        nst = [l.strip() for l in head if re.match(r"^\s*nstates\s*=", l)]
        cls = sorted(set(re.findall(r"from (\w+) import (\w+)", txt)))
        xc = sorted(set(re.findall(r"XC functional[s]? =\s*(\S+)", txt)))
        gene = sorted(set(re.findall(r"\b(ALDA0|Multicollinear)\b", txt)))
        n_ee = len(re.findall(r"^excited energy", txt, re.M))
        labels = re.findall(r"^(ALDA0_\w+|Multicollinear_\w+)(.*)$", txt, re.M)
        mols = sorted(set(re.findall(r"^\[INPUT\] num\. atoms\s*=\s*(\d+)", txt, re.M)))
        print(f"\n{os.path.basename(path)}  ({os.path.getsize(path)} B, {len(lines)} lines)")
        print(f"   class={cls}  funs={funs}  {nst}")
        print(f"   XC={xc}   'excited energy' lines={n_ee}")
        lab = [l for l, tail in labels]
        print(f"   labels: {len(lab)}  first={lab[:2]}  counts={ {l: lab.count(l) for l in sorted(set(lab))} }")
        print(f"   mol inline in label: {sum(1 for _, t in labels if t.strip())}")
        print(f"   distinct num-atoms: {len(mols)}")
