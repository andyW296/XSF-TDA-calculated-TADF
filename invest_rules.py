"""Which state-2/3 rule reproduces the Peng USF-TDA columns?"""
import pandas as pd
import mr_tadf_lib as L
import invest_audit as IA

SF, SFSTATE, df = IA.SF, IA.SFSTATE, IA.df
row_of = {}
for r in range(3, df.shape[0]):
    v = df.iat[r, 0]
    if pd.notna(v):
        row_of[int(v)] = r

RULES = {
    "smaller D<S^2> of 2,3": lambda d: L.pick_state_2_3(d, keep="smaller")[0],
    "larger  D<S^2> of 2,3": lambda d: L.pick_state_2_3(d, keep="larger")[0],
    "state 2 always": lambda d: next((e for n, e, s in d if n == 2), None),
    "state 3 always": lambda d: next((e for n, e, s in d if n == 3), None),
    "1st state from 2 with D<S^2><-0.3":
        lambda d: next((e for n, e, s in d if n >= 2 and s is not None and s < -0.3), None),
}
print(f"{'col':>4} {'method/var/fun':<34} " + " ".join(f"{k[:22]:>23}" for k in RULES))
for cidx, (kind, meth, variant, fun) in L.COLUMNS.items():
    if kind != "SF" or meth != "USF-TDA":
        continue
    letter = L.letter(cidx)
    d = SF.get((meth, variant, fun), {})
    det = SFSTATE.get((meth, variant, fun), {})
    counts = {k: [0, 0] for k in RULES}          # [ok, compared]
    for m, rr in d.items():
        if m not in row_of:
            continue
        ex = df.iat[row_of[m], cidx]
        if pd.isna(ex):
            continue
        detail = det[m][0][0] if m in det else []
        for name, fn in RULES.items():
            v = fn(detail)
            if v is not None:
                counts[name][1] += 1
                if abs(v - float(ex)) <= 5e-4:
                    counts[name][0] += 1
    print(f"{letter:>4} {meth+'/'+variant+'/'+fun:<34} " +
          " ".join(f"{counts[k][0]:>10}/{counts[k][1]:<12}" for k in RULES))
