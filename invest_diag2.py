"""Dump the spin-flip spectra behind the disagreeing Peng cells."""
import pandas as pd
import mr_tadf_lib as L
import invest_audit as IA

TD, SF, SFSTATE, USTD, df = IA.TD, IA.SF, IA.SFSTATE, IA.USTD, IA.df
row_of = {}
for r in range(3, df.shape[0]):
    v = df.iat[r, 0]
    if pd.notna(v):
        row_of[int(v)] = r

CASES = [("N", 22), ("N", 31), ("P", 36), ("V", 6), ("V", 30), ("V", 34), ("V", 37),
         ("V", 39), ("V", 44), ("AD", 6), ("AE", 37), ("O", 44)]
letters = {L.letter(c): (c, v) for c, v in L.COLUMNS.items()}

for letter, mol in CASES:
    cidx, (kind, meth, variant, fun) = letters[letter]
    ex = df.iat[row_of[mol], cidx]
    print("=" * 104)
    print(f"col {letter} = {meth} {variant} {fun} | molecule {mol} | excel = {ex}")
    d = SF.get((meth, variant, fun), {}).get(mol, [])
    det = SFSTATE.get((meth, variant, fun), {}).get(mol, [])
    for k, (nums, src) in enumerate(d):
        print(f"   {src:<26} " + " ".join(f"{x:>8.4f}" for x in nums))
        if k < len(det):
            dd = det[k][0]
            print(f"   {'':<26} " + " ".join(f"{v:>8.4f}" for _, v, _ in dd) + "  <- per-state E")
            print(f"   {'':<26} " + " ".join(
                (f"{s:>8.2f}" if s is not None else f"{'--':>8}") for _, _, s in dd)
                + "  <- D<S^2>")
    print(f"   excel value {ex} matches: ", end="")
    hits = []
    for nums, src in d:
        for i, x in enumerate(nums):
            if abs(x - float(ex)) <= 5e-4:
                hits.append(f"{src} state {i+1}")
    print(", ".join(hits) if hits else "nothing in these arrays")
