"""Print the BHHLYP USF-TD-DFT states (2nd/3rd + <S^2>) of all 35 molecules."""
import mr_tadf_lib as L

TD, SF, SFSTATE, USTD, inv = L.collect()
df, names = L.load_excel()

lines = []
hdr = (f"{'#':>3} {'molecule':<22} {'state1 E/eV':>11} {'<S^2>':>6} "
       f"{'state2 E/eV':>11} {'<S^2>':>6} {'state3 E/eV':>11} {'<S^2>':>6} "
       f"{'keep':>4} {'rule value':>10} {'excel U':>9} {'diff':>8}")
lines.append("BHHLYP  USF-TD-DFT : excited states 2 and 3 and their <S^2>")
lines.append("rule: of states 2 and 3 keep the one with the SMALLER <S^2>")
lines.append("source: alda0_mc_bhhlyp.log (= bhhlyp_1_10.log) for 01-10, mc_bhhlyp.log for 11-35")
lines.append("")
lines.append(hdr)
lines.append("-" * len(hdr))

n_ok = 0
for m in range(1, 36):
    ex = df.iat[3 + m - 1, 20]
    rr = USTD.get("BHHLYP", {}).get(m)
    if not rr:
        lines.append(f"{m:>3} {str(names[m-1]):<22}   -- no USF-TD-DFT data --")
        continue
    states, src = rr[0]
    d = {n: (e, s) for n, e, s in states}
    v, n, spin = L.pick_state_2_3(states, keep="smaller")

    def f(k):
        return f"{d[k][0]:>11.5f} {d[k][1]:>6.3f}" if k in d else f"{'--':>11} {'--':>6}"
    diff = v - float(ex) if (v is not None) else float("nan")
    ok = "" if abs(diff) <= 5e-4 else "  <-- differs"
    if abs(diff) <= 5e-4:
        n_ok += 1
    lines.append(f"{m:>3} {str(names[m-1]):<22} {f(1)} {f(2)} {f(3)} {n:>4} "
                 f"{v:>10.5f} {float(ex):>9.5f} {diff:>+8.5f}{ok}")

lines.append("-" * len(hdr))
lines.append(f"molecules listed: 35, cells reproduced by the rule: {n_ok}")
out = "\n".join(lines)
open("usf_tddft_bhhlyp_states.txt", "w", encoding="utf-8").write(out + "\n")
print(out)
