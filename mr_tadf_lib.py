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
# E (4), M (12) and AC (28) = USF-TD-DFT for B3LYP / PBE0 / CAM-B3LYP are excluded:
# those calculations do not exist in this repository (user instruction).
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
    20: ("USTD", "USF-TD-DFT", None, "BHHLYP"),
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
EXCLUDED = {4: "E", 12: "M", 28: "AC"}      # USF-TD-DFT / B3LYP, PBE0, CAM-B3LYP
NO_SOURCE = {34: "CV1", 35: "CV1", 36: "CV1", 37: "CV1", 38: "CV1", 39: "CV1",
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
# 'Excited State   3:      0.26345 eV   <S^2>: 0.2223'  (spin-flip TD-DFT blocks)
RE_USTATE = re.compile(r"^Excited State\s+(\d+):\s+(-?\d+\.\d+)\s+eV\s+<S\^2>:\s+(\d+\.\d+)")
RE_USLABEL = re.compile(r"^USF[- ]TD-DFT\b")
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
    """Sections per gto.M call: geometry, functional, TD analyses, SF blocks.

    The spin-flip blocks are labelled with the method they actually belong to, decided
    from their own output (SA_SF_TDA prints D<S^2> per state, XSF_TDA prints CV(1): ...
    instead), because the logs are not always stored in the folder of their method.
    """
    sections, cur, pending, variant = [], None, None, None
    buf = None
    versions = []
    text = open(path, encoding="utf-8", errors="ignore").read()
    has_sf = "SA_SF_TDA" in text or "SF_XTDA" in text
    has_xsf = "XSF_TDA" in text
    file_default = "USF-TDA" if (has_sf and not has_xsf) else (
        "XSF-TDA" if (has_xsf and not has_sf) else "XSF-TDA")
    # some runs print no ALDA0 / Multicollinear label at all; the file name then says it
    low = os.path.basename(path).lower()
    hint = "MC" if ("mc" in low and "alda0" not in low) else (
        "ALDA0" if ("alda0" in low and "mc" not in low) else None)
    for raw in text.splitlines():
        s = raw.strip()
        m = RE_PYSCF.match(s)
        if m:
            versions.append(m.group(1))

        m = RE_NUMATOMS.match(s)
        if m:
            if cur is not None:
                _fix_variants(cur, hint)
            cur = {"natoms": int(m.group(1)), "atoms": [], "xc": None, "td": [],
                   "sf": [], "sdetail": [], "us": []}
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
        if m:
            # a section may run the SCF with one functional and then switch mf.xc before the
            # spin-flip / TD step, so the LAST 'XC functionals' line is the operative one
            cur["xc"] = m.group(1).lower()
            continue

        # ---- 'excited energy [ ... ]' array (may wrap over lines) ------------
        if buf is not None:
            buf += " " + s
            if "]" in s:
                nums = [float(x) for x in re.findall(r"-?\d+\.\d+", buf)]
                cur["sf"].append({"variant": variant, "nums": nums, "detail": [],
                                  "usf": False, "xsf": False})
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
                cur["sf"].append({"variant": variant, "nums": nums, "detail": [],
                                  "usf": False, "xsf": False})
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
        # ---- per-state detail lines of an SF analyze() block -----------------
        if s.startswith("CV(") and cur["sf"]:
            cur["sf"][-1]["xsf"] = True          # XSF_TDA prints CV(1): ... instead
            continue
        m = RE_ES.match(s)
        if m:
            n = int(m.group(1))
            d = RE_DS2.search(s)
            sv = float(d.group(1)) if d else None
            # a new block starts when no block is open, or when the state counter
            # restarts at 1 (the log then prints only the state table, no array)
            if not cur["sf"] or (n == 1 and cur["sf"][-1]["detail"]):
                cur["sf"].append({"variant": variant, "nums": None, "detail": [],
                                  "usf": False, "xsf": False})
            cur["sf"][-1]["detail"].append((n, float(m.group(2)), sv))
            if sv is not None:
                cur["sf"][-1]["usf"] = True
            continue
        # ---- spin-flip TD-DFT block: 'Excited State   N:  E eV  <S^2>: x' ----
        if RE_USLABEL.match(s):
            cur["us"].append([])
            continue
        m = RE_USTATE.match(s)
        if m:
            if not cur["us"]:
                cur["us"].append([])
            cur["us"][-1].append((int(m.group(1)), float(m.group(2)), float(m.group(3))))
    if cur is not None:
        _fix_variants(cur, hint)
    return sections, versions, file_default


def _fix_variants(sec, hint=None):
    """Assign a variant to blocks whose label was not printed: the code runs ALDA0
    first and the multicollinear variant second inside one gto.M section, unless the
    file name says the whole run is a single variant (hint).  Also fill in the
    eigenvalue list from the per-state table when no 'excited energy' array was printed.
    """
    used = [b["variant"] for b in sec["sf"] if b["variant"]]
    # some scripts print 'ALDA0' in front of both blocks (copy-paste bug); the code always
    # runs ALDA0 first and the multicollinear variant second, so fix that pair by position
    if (len(sec["sf"]) == 2 and sec["sf"][0]["variant"] == "ALDA0"
            and sec["sf"][1]["variant"] == "ALDA0"):
        sec["sf"][0]["variant"] = "ALDA0"
        sec["sf"][1]["variant"] = "MC"
        used = ["ALDA0", "MC"]
    if hint and not used:                      # unlabelled single-variant run
        for b in sec["sf"]:
            b["variant"] = hint
        used = [hint]
    for b in sec["sf"]:
        if b["variant"] is None:
            if "ALDA0" not in used:
                b["variant"] = "ALDA0"
            elif "MC" not in used:
                b["variant"] = "MC"
            else:
                b["variant"] = "EXTRA"
            used.append(b["variant"])
        if b["nums"] is None:
            b["nums"] = [e for _, e, _ in sorted(b["detail"])]


def pick_state_2_3(detail, keep="smaller"):
    """Comparison rule: look at excited states 2 and 3 only and keep the one with
    either the smaller or the larger spin expectation value.

    detail = [(state_number, energy_eV, spin_value or None), ...]
      USF-TDA      -> keep="smaller"  on D<S^2>
      USF-TD-DFT   -> keep="larger"   on <S^2>
    Returns (value, state_number, spin_value) or (None, None, None) when the spin
    expectation values are not printed.
    """
    cand = {n: (e, s) for n, e, s in detail if n in (2, 3) and s is not None}
    if not cand:
        return None, None, None
    pick = min if keep == "smaller" else max
    n = pick(cand, key=lambda k: cand[k][1])
    return cand[n][0], n, cand[n][1]


def load_excel():
    df = pd.ExcelFile(XLSX).parse("MR-TADF", header=None)
    names = df.iloc[3:, 0].tolist()
    return df, names


def collect(folders=None, geo2mol=None):
    """Return (TD, SF, SFSTATE, USTD, inventory) for the given log folders."""
    if folders is None:
        folders = ["MR-TADF/TDDFT", "MR-TADF/USF-TDA", "MR-TADF/XSF-TDA"]
    if geo2mol is None:
        geo2mol = build_geo_map()
    TD = collections.defaultdict(lambda: collections.defaultdict(list))
    SF = collections.defaultdict(lambda: collections.defaultdict(list))
    SFSTATE = collections.defaultdict(lambda: collections.defaultdict(list))
    USTD = collections.defaultdict(lambda: collections.defaultdict(list))
    inventory = []
    for folder in folders:
        tag = os.path.basename(folder)
        for path in sorted(glob.glob(os.path.join(folder, "*"))):
            if not os.path.isfile(path):
                continue
            name = os.path.basename(path)
            sections, versions, file_default = parse_log(path)
            last_xc, n_id, n_sf, n_td, n_us = None, 0, 0, 0, 0
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
                for k, blk in enumerate(sec["sf"]):
                    # the method is decided from the block's own output
                    if blk["usf"]:
                        blk_meth = "USF-TDA"
                    elif blk["xsf"]:
                        blk_meth = "XSF-TDA"
                    else:
                        blk_meth = file_default
                    SF[(blk_meth, blk["variant"], fun)][mol].append((blk["nums"], name))
                    SFSTATE[(blk_meth, blk["variant"], fun)][mol].append(
                        (blk["detail"], name, k))
                    per[(blk_meth + "/" + blk["variant"], fun)] += 1
                    n_sf += 1
                for states in sec["us"]:
                    USTD[fun][mol].append((states, name))
                    per[("USF-TDDFT", fun)] += 1
                    n_us += 1
                if len(sec["sf"]) == 1:
                    weird.append((mol, sec["sf"][0]["variant"]))
                if not sec["td"] and not sec["sf"] and not sec["us"]:
                    per[("none", fun)] += 1
            inventory.append({"folder": tag, "log": name, "sections": len(sections),
                              "identified": n_id, "td": n_td, "sf": n_sf, "us": n_us,
                              "versions": sorted(set(versions)), "per": dict(per),
                              "single_array": weird})
    return TD, SF, SFSTATE, USTD, inventory
