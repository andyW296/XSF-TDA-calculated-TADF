"""Shared parsing library for the MR-TADF log folders."""
import os
import re
import glob
import collections
import pandas as pd

XYZ_DIR = "MR-TADF/xyz"
XLSX = "TADF_all.xlsx"
FUN_NORM = {"b3lyp": "B3LYP", "pbe0": "PBE0", "bhandhlyp": "BHHLYP",
            "cam-b3lyp": "CAM-B3LYP", "cam_b3lyp": "CAM-B3LYP", "camb3lyp": "CAM-B3LYP"}

# excel column index -> (kind, folder-or-method, variant, functional)
COLUMNS = {
    2:  ("TD", "TDA",     None,   "B3LYP"),
    3:  ("TD", "TD-DFT",  None,   "B3LYP"),
    5:  ("SF", "USF-TDA", "ALDA0", "B3LYP"),
    6:  ("SF", "XSF-TDA", "ALDA0", "B3LYP"),
    7:  ("SF", "USF-TDA", "MC",    "B3LYP"),
    8:  ("SF", "XSF-TDA", "MC",    "B3LYP"),
    10: ("TD", "TDA",     None,   "PBE0"),
    11: ("TD", "TD-DFT",  None,   "PBE0"),
    13: ("SF", "USF-TDA", "ALDA0", "PBE0"),
    14: ("SF", "XSF-TDA", "ALDA0", "PBE0"),
    15: ("SF", "USF-TDA", "MC",    "PBE0"),
    16: ("SF", "XSF-TDA", "MC",    "PBE0"),
    18: ("TD", "TDA",     None,   "BHHLYP"),
    19: ("TD", "TD-DFT",  None,   "BHHLYP"),
    21: ("SF", "USF-TDA", "ALDA0", "BHHLYP"),
    22: ("SF", "XSF-TDA", "ALDA0", "BHHLYP"),
    23: ("SF", "USF-TDA", "MC",    "BHHLYP"),
    24: ("SF", "XSF-TDA", "MC",    "BHHLYP"),
    26: ("TD", "TDA",     None,   "CAM-B3LYP"),
    27: ("TD", "TD-DFT",  None,   "CAM-B3LYP"),
    29: ("SF", "USF-TDA", "ALDA0", "CAM-B3LYP"),
    30: ("SF", "XSF-TDA", "ALDA0", "CAM-B3LYP"),
    31: ("SF", "USF-TDA", "MC",    "CAM-B3LYP"),
    32: ("SF", "XSF-TDA", "MC",    "CAM-B3LYP"),
}
NO_SOURCE = {4: "USF-TD-DFT/B3LYP", 12: "USF-TD-DFT/PBE0", 20: "USF-TD-DFT/BHHLYP",
             28: "USF-TD-DFT/CAM-B3LYP",
             34: "CV1", 35: "CV1", 36: "CV1", 37: "CV1", 38: "CV1", 39: "CV1",
             40: "CV1", 41: "CV1",
             43: "DS2", 44: "DS2", 45: "DS2", 46: "DS2", 47: "DS2", 48: "DS2",
             49: "DS2", 50: "DS2", 51: "DS2", 52: "DS2", 53: "DS2", 54: "DS2",
             56: "OO1"}


def letter(i):
    return chr(65 + i) if i < 26 else "A" + chr(65 + i - 26)


RE_NUMATOMS = re.compile(r"^\[INPUT\] num\. atoms\s*=\s*(\d+)")
RE_INPUT_ATOM = re.compile(r"^\[INPUT\]\s*(\d+)\s+([A-Za-z]{1,2})\s+"
                           r"(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+AA\b")
RE_XC = re.compile(r"XC functionals\s*=\s*(\S+)")
RE_BANNER = re.compile(r"<class 'pyscf\.tdscf\.rks\.(TDA|TDDFT)'>")
RE_SPIN = re.compile(r"\*\* (Singlet|Triplet) excitation energies")
RE_ES1 = re.compile(r"^Excited State\s+1:\s+(-?\d+\.\d+)\s+eV")
RE_ES = re.compile(r"^Excited state\s+(\d+)\s+(-?\d+\.\d+)\s+eV")
RE_DS2 = re.compile(r"D<S\^2>=(-?\d+\.\d+)")
RE_EXEN = re.compile(r"^excited energy\s+\[")
RE_PYSCF = re.compile(r"^PySCF version\s+(\S+)")
# substring match: the labels are printed as 'ALDA0_b3lyp' / 'Multicollinear_b3lyp',
# so a trailing \b would never match before the underscore
RE_HAS_ALDA0 = re.compile(r"ALDA0")
RE_HAS_MC = re.compile(r"Multicollinear")
RE_LEGEND = re.compile(r"method=0.*ALDA0.*method=1", re.I)


def read_xyz(path):
    atoms = []
    for line in open(path, encoding="utf-8", errors="ignore").read().splitlines()[2:]:
        p = line.split()
        if len(p) >= 4:
            atoms.append((p[0], round(float(p[1]), 2), round(float(p[2]), 2),
                          round(float(p[3]), 2)))
    return tuple(sorted(atoms))


def build_geo_map():
    geo2mol = {}
    for d in sorted(os.listdir(XYZ_DIR)):
        p = os.path.join(XYZ_DIR, d, "geom.xyz")
        if os.path.isfile(p):
            geo2mol[read_xyz(p)] = int(d)
    return geo2mol


def parse_log(path):
    """Sections per gto.M call: geometry, functional, TD analyses, SF arrays."""
    sections, cur, pending, variant = [], None, None, None
    buf = None
    versions = []
    for raw in open(path, encoding="utf-8", errors="ignore"):
        s = raw.strip()
        m = RE_PYSCF.match(s)
        if m:
            versions.append(m.group(1))

        m = RE_NUMATOMS.match(s)
        if m:
            if cur is not None:
                _fix_variants(cur)
            cur = {"natoms": int(m.group(1)), "atoms": [], "xc": None, "td": [],
                   "sf": [], "sdetail": []}
            sections.append(cur)
            pending, variant, buf = None, None, None
            continue
        if cur is None:
            continue

        m = RE_INPUT_ATOM.match(s)
        if m:
            cur["atoms"].append((m.group(2), round(float(m.group(3)), 2),
                                 round(float(m.group(4)), 2), round(float(m.group(5)), 2)))
            continue
        m = RE_XC.search(s)
        if m and cur["xc"] is None:
            cur["xc"] = m.group(1).lower()
            continue

        # ---- 'excited energy [ ... ]' array (may wrap over lines) ------------
        if buf is not None:
            buf += " " + s
            if "]" in s:
                nums = [float(x) for x in re.findall(r"-?\d+\.\d+", buf)]
                cur["sf"].append([variant, nums])
                cur["sdetail"].append([])
                buf = None
            continue
        if not RE_LEGEND.search(s):
            if RE_HAS_MC.search(s):
                variant = "MC"
            elif RE_HAS_ALDA0.search(s):
                variant = "ALDA0"
        if RE_EXEN.match(s):
            if "]" in s:
                nums = [float(x) for x in re.findall(r"-?\d+\.\d+", s)]
                cur["sf"].append([variant, nums])
                cur["sdetail"].append([])
            else:
                buf = s
            continue

        # ---- TD-DFT folder output -------------------------------------------
        m = RE_BANNER.search(s)
        if m:
            pending = {"method": m.group(1), "spin": None}
            continue
        if pending is not None:
            m = RE_SPIN.search(s)
            if m and pending["spin"] is None:
                pending["spin"] = m.group(1)
                continue
            m = RE_ES1.match(s)
            if m and pending["spin"] is not None:
                cur["td"].append((pending["method"], pending["spin"], float(m.group(1))))
                pending = None
                continue
        # ---- per-state detail lines (shift relative to the SF reference) -----
        m = RE_ES.match(s)
        if m:
            d = RE_DS2.search(s)
            if cur["sdetail"]:
                cur["sdetail"][-1].append((int(m.group(1)), float(m.group(2)),
                                           float(d.group(1)) if d else None))
    if cur is not None:
        _fix_variants(cur)
    return sections, versions


def _fix_variants(sec):
    """Assign a variant to arrays whose label was not printed: the code always runs
    ALDA0 first and the multicollinear variant second inside one gto.M section."""
    used = [a[0] for a in sec["sf"] if a[0]]
    for a in sec["sf"]:
        if a[0] is None:
            if "ALDA0" not in used:
                a[0] = "ALDA0"
            elif "MC" not in used:
                a[0] = "MC"
            else:
                a[0] = "EXTRA"
            used.append(a[0])


def usf_first_singlet(nums, detail):
    """USF-TDA: energy of the lowest non-S0 state with singlet character.

    The analyze output prints D<S^2> per state; a clean triplet has D<S^2> ~ 0 and a
    singlet (including the open-shell/contaminated ones) has D<S^2> < -0.3.
    Returns (value, state_index0, pure) or (None, None, False) if unknowable."""
    if detail:
        for n, e, ds in detail:
            if n <= 1 or ds is None:
                continue
            if ds < -0.3:
                return e, n - 1, True
    return None, None, False


def load_excel():
    df = pd.ExcelFile(XLSX).parse("MR-TADF", header=None)
    names = df.iloc[3:, 0].tolist()
    return df, names


def collect():
    """Return (TD, SF, SFSTATE, inventory)."""
    geo2mol = build_geo_map()
    TD = collections.defaultdict(lambda: collections.defaultdict(list))
    SF = collections.defaultdict(lambda: collections.defaultdict(list))
    SFSTATE = collections.defaultdict(lambda: collections.defaultdict(list))
    inventory = []
    for folder in ["MR-TADF/TDDFT", "MR-TADF/USF-TDA", "MR-TADF/XSF-TDA"]:
        tag = os.path.basename(folder)
        for path in sorted(glob.glob(os.path.join(folder, "*"))):
            if not os.path.isfile(path):
                continue
            name = os.path.basename(path)
            sections, versions = parse_log(path)
            last_xc, n_id, n_sf, n_td = None, 0, 0, 0
            per = collections.Counter()
            weird = []
            for sec in sections:
                if sec["xc"] is None:
                    sec["xc"] = last_xc
                else:
                    last_xc = sec["xc"]
                mol = geo2mol.get(tuple(sorted(sec["atoms"])))
                if mol is None:
                    continue
                n_id += 1
                fun = FUN_NORM.get(sec["xc"] or "", sec["xc"] or "?")
                for meth, spin, e in sec["td"]:
                    m = {"TDDFT": "TD-DFT", "TDA": "TDA"}.get(meth, meth)
                    TD[(m, fun, spin)][mol].append((e, name))
                    n_td += 1
                for k, (variant, nums) in enumerate(sec["sf"]):
                    SF[(tag, variant, fun)][mol].append((nums, name))
                    det = sec["sdetail"][k] if k < len(sec["sdetail"]) else []
                    SFSTATE[(tag, variant, fun)][mol].append((det, name, k))
                    per[(variant, fun)] += 1
                    n_sf += 1
                if len(sec["sf"]) == 1:
                    weird.append((mol, sec["sf"][0][0]))
                if not sec["td"] and not sec["sf"]:
                    per[("none", fun)] += 1
            inventory.append({"folder": tag, "log": name, "sections": len(sections),
                              "identified": n_id, "td": n_td, "sf": n_sf,
                              "versions": sorted(set(versions)), "per": dict(per),
                              "single_array": weird})
    return TD, SF, SFSTATE, inventory
