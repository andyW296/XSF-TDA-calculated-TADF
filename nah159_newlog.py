import re
import collections

p = "NAH159/XSF-TDA/tddft_bhhlyp.log2"
t = open(p, encoding="utf-8", errors="ignore").read()
lines = t.splitlines()
c = collections.Counter()
method = None
for l in lines:
    s = l.strip()
    m = re.search(r"<class 'pyscf\.tdscf\.rks\.(TDA|TDDFT)'>", s)
    if m:
        method = m.group(1)
    m = re.match(r"\*\* (Singlet|Triplet) excitation energies", s)
    if m and method:
        c[(method, m.group(1))] += 1
print(p, dict(c))
print("funs:", sorted(set(re.findall(r"funs? = (.+)", t)))[:3])
print("mols slice:", [l.strip() for l in lines[:200] if "mols" in l][:3])
